#!/usr/bin/env python3
"""
Colorado Snowfall Dashboard - Data Fetcher

Fetches SNOTEL data from NRCS AWDB REST API, computes statewide and per-basin
aggregations, historical percentile envelope, and catch-up analysis.

Usage:
    pip install requests
    python fetch_data.py

Output: data.js (JavaScript file for the dashboard)
"""

import requests
import json
import datetime
import time
import sys
from collections import defaultdict

# --- Configuration ---
BASE_URL = "https://wcc.sc.egov.usda.gov/awdbRestApi/services/v1"
BATCH_SIZE = 25  # stations per API call
RATE_LIMIT_SLEEP = 0.3  # seconds between API calls
HISTORICAL_YEARS = 35  # WY 1991 through present (SNOTEL data reliable from ~1987)

# Water year: Oct 1 - Sep 30. WY 2026 started Oct 1, 2025.
TODAY = datetime.date.today()
if TODAY.month >= 10:
    CURRENT_WY = TODAY.year + 1
else:
    CURRENT_WY = TODAY.year

WY_START = datetime.date(CURRENT_WY - 1, 10, 1)

# Colorado basin mapping by HUC4 prefix
BASIN_MAP = {
    "1019": "South Platte",
    "1018": "North Platte / Laramie",
    "1401": "Upper Colorado",
    "1402": "Gunnison",
    "1403": "Dolores / San Miguel",
    "1405": "Yampa / White",
    "1408": "San Juan / Animas",
    "1102": "Arkansas",
    "1301": "Rio Grande",
}


def api_get(endpoint, params, retries=3):
    """Make a GET request to the AWDB API with retries."""
    url = f"{BASE_URL}/{endpoint}"
    for attempt in range(retries):
        try:
            resp = requests.get(url, params=params, timeout=60)
            if resp.status_code == 200:
                return resp.json()
            elif resp.status_code in (429, 503):
                wait = (attempt + 1) * 5
                print(f"  Rate limited ({resp.status_code}), waiting {wait}s...")
                time.sleep(wait)
            else:
                print(f"  API error {resp.status_code}: {resp.text[:200]}")
                if attempt < retries - 1:
                    time.sleep(2)
        except requests.exceptions.RequestException as e:
            print(f"  Request failed: {e}")
            if attempt < retries - 1:
                time.sleep(2)
    return None


def get_stations():
    """Fetch all active Colorado SNOTEL stations."""
    print("Fetching CO SNOTEL stations...")
    # API ignores filter params, so we fetch all and filter client-side
    data = api_get("stations", {})
    if not data:
        print("ERROR: Could not fetch stations")
        sys.exit(1)

    stations = []
    for s in data:
        triplet = s.get("stationTriplet", "")
        # Filter to Colorado SNOTEL stations only
        if ":CO:SNTL" not in triplet:
            continue
        # Skip inactive stations (endDate in the past)
        end_date = s.get("endDate", "")
        if end_date and "2100" not in end_date:
            continue
        huc = s.get("huc", "")
        basin = "Other"
        for prefix, name in BASIN_MAP.items():
            if huc.startswith(prefix):
                basin = name
                break
        stations.append({
            "triplet": triplet,
            "name": s.get("name", ""),
            "elevation": s.get("elevation", 0),
            "lat": s.get("latitude", 0),
            "lon": s.get("longitude", 0),
            "huc": huc,
            "basin": basin,
        })

    print(f"  Found {len(stations)} active CO SNOTEL stations")
    return stations


def batch_station_triplets(stations):
    """Yield batches of station triplets."""
    triplets = [s["triplet"] for s in stations]
    for i in range(0, len(triplets), BATCH_SIZE):
        yield triplets[i:i + BATCH_SIZE]


def fetch_swe_batch(triplet_batch, begin_date, end_date, with_median=True):
    """Fetch daily SWE data for a batch of stations."""
    params = {
        "stationTriplets": ",".join(triplet_batch),
        "elements": "WTEQ",
        "duration": "DAILY",
        "beginDate": begin_date,
        "endDate": end_date,
        "periodRef": "START",
        "getFlags": "false",
    }
    if with_median:
        params["centralTendencyType"] = "MEDIAN"
    return api_get("data", params)


def fetch_all_stations(stations, begin_date, end_date, with_median=True, label=""):
    """Fetch SWE data for all stations in batches."""
    all_data = []
    batches = list(batch_station_triplets(stations))
    for i, batch in enumerate(batches):
        if label:
            print(f"  {label} batch {i+1}/{len(batches)}...")
        result = fetch_swe_batch(batch, begin_date, end_date, with_median)
        if result:
            all_data.extend(result)
        time.sleep(RATE_LIMIT_SLEEP)
    return all_data


def aggregate_by_date(all_station_data, station_filter=None):
    """
    Compute mean SWE and mean median by date across stations.
    Returns dict: date_str -> {swe, median, count}
    """
    date_vals = defaultdict(lambda: {"swe_vals": [], "med_vals": []})

    for station_result in all_station_data:
        if not station_result:
            continue
        triplet = station_result.get("stationTriplet", "")
        if station_filter and triplet not in station_filter:
            continue
        for data_block in station_result.get("data", []):
            for v in data_block.get("values", []):
                date_str = v.get("date", "")
                # Handle date formats like "2025-10-01 00:00" -> "2025-10-01"
                if " " in date_str:
                    date_str = date_str.split(" ")[0]
                val = v.get("value")
                med = v.get("median")
                if val is not None:
                    date_vals[date_str]["swe_vals"].append(val)
                if med is not None:
                    date_vals[date_str]["med_vals"].append(med)

    result = {}
    for date_str in sorted(date_vals.keys()):
        d = date_vals[date_str]
        swe_vals = d["swe_vals"]
        med_vals = d["med_vals"]
        if len(swe_vals) >= 10:  # require at least 10 stations reporting
            result[date_str] = {
                "swe": round(sum(swe_vals) / len(swe_vals), 2),
                "median": round(sum(med_vals) / len(med_vals), 2) if med_vals else 0,
                "count": len(swe_vals),
            }
    return result


def date_to_doy(date_str, wy_start_year):
    """Convert date string to day-of-water-year (Oct 1 = day 1)."""
    d = datetime.date.fromisoformat(date_str)
    wy_start = datetime.date(wy_start_year, 10, 1)
    return (d - wy_start).days + 1


def doy_to_date_str(doy, wy_start_year):
    """Convert day-of-water-year back to date string."""
    wy_start = datetime.date(wy_start_year, 10, 1)
    d = wy_start + datetime.timedelta(days=doy - 1)
    return d.isoformat()


def compute_envelope(historical_curves, max_doy):
    """Compute percentile bands from historical year curves."""
    envelope = {
        "min": [], "p10": [], "p25": [], "p75": [], "p90": [], "max": []
    }
    for doy in range(1, max_doy + 1):
        vals = []
        for yr_curve in historical_curves.values():
            if doy in yr_curve:
                vals.append(yr_curve[doy])
        if len(vals) >= 5:
            vals.sort()
            n = len(vals)
            envelope["min"].append(round(vals[0], 2))
            envelope["max"].append(round(vals[-1], 2))
            envelope["p10"].append(round(vals[max(0, int(n * 0.1))], 2))
            envelope["p25"].append(round(vals[max(0, int(n * 0.25))], 2))
            envelope["p75"].append(round(vals[min(n - 1, int(n * 0.75))], 2))
            envelope["p90"].append(round(vals[min(n - 1, int(n * 0.9))], 2))
        else:
            for key in envelope:
                envelope[key].append(None)
    return envelope


def main():
    print("=" * 60)
    print("  Colorado Snowfall Dashboard - Data Fetcher")
    print(f"  Water Year {CURRENT_WY} | Today: {TODAY}")
    print("=" * 60)

    # -------------------------------------------------------
    # Step 1: Get stations
    # -------------------------------------------------------
    stations = get_stations()
    station_triplets = set(s["triplet"] for s in stations)

    # Build basin -> triplet mapping
    basin_triplets = defaultdict(set)
    for s in stations:
        basin_triplets[s["basin"]].add(s["triplet"])

    # -------------------------------------------------------
    # Step 2: Fetch current water year SWE data with medians
    # -------------------------------------------------------
    print(f"\nFetching current season data ({WY_START} to {TODAY})...")
    current_data = fetch_all_stations(
        stations, WY_START.isoformat(), TODAY.isoformat(),
        with_median=True, label="Current season"
    )
    current_agg = aggregate_by_date(current_data)
    print(f"  Got {len(current_agg)} days of statewide data")

    # -------------------------------------------------------
    # Step 3: Fetch historical water years
    # (We derive the full-season median curve from this data)
    # -------------------------------------------------------
    hist_start_wy = CURRENT_WY - HISTORICAL_YEARS
    print(f"\nFetching {HISTORICAL_YEARS} historical water years ({hist_start_wy}-{CURRENT_WY - 1})...")

    historical_curves = {}  # {wy_year: {doy: mean_swe}}
    for wy in range(hist_start_wy, CURRENT_WY):
        wy_begin = f"{wy - 1}-10-01"
        wy_end = f"{wy}-06-30"
        print(f"  WY {wy}...", end="", flush=True)

        all_hist = fetch_all_stations(
            stations, wy_begin, wy_end,
            with_median=False, label=""
        )
        hist_agg = aggregate_by_date(all_hist)

        curve = {}
        for date_str, vals in hist_agg.items():
            doy = date_to_doy(date_str, wy - 1)
            curve[doy] = vals["swe"]
        historical_curves[str(wy)] = curve
        print(f" {len(curve)} days")

    # -------------------------------------------------------
    # Step 5: Compute historical envelope AND full-season median curve
    # -------------------------------------------------------
    print("\nComputing historical envelope and median curve...")
    max_doy = 274  # Oct 1 through Jun 30

    envelope = compute_envelope(historical_curves, max_doy)

    # Build the full-season median curve from historical data (median of yearly means)
    # This approximates the 30-year NRCS median well enough for visualization
    full_median_curve = {}
    date_labels = []
    for doy in range(1, max_doy + 1):
        d = WY_START + datetime.timedelta(days=doy - 1)
        date_str = d.isoformat()
        date_labels.append(date_str)
        vals = [historical_curves[yr][doy] for yr in historical_curves if doy in historical_curves[yr]]
        if len(vals) >= 5:
            vals.sort()
            # Use true median (middle value)
            n = len(vals)
            if n % 2 == 0:
                med = (vals[n // 2 - 1] + vals[n // 2]) / 2
            else:
                med = vals[n // 2]
            full_median_curve[date_str] = round(med, 2)

    # For dates we have actual API median data (current season), use those instead
    # as they're the official 30-year 1991-2020 medians
    for date_str, agg in current_agg.items():
        if agg["median"] > 0:
            full_median_curve[date_str] = agg["median"]

    print(f"  Full-season median curve: {len(full_median_curve)} days")

    # -------------------------------------------------------
    # Step 6: Compute summary stats and catch-up analysis
    # -------------------------------------------------------
    print("Computing catch-up analysis...")

    sorted_dates = sorted(current_agg.keys())
    if not sorted_dates:
        print("ERROR: No current season data")
        sys.exit(1)

    latest_date = sorted_dates[-1]
    current_swe = current_agg[latest_date]["swe"]
    median_today = current_agg[latest_date]["median"]
    deficit = round(median_today - current_swe, 2)
    pct_of_median = round((current_swe / median_today) * 100, 1) if median_today > 0 else 0

    # Find median peak from full-season median curve
    median_peak_swe = 0
    median_peak_date = None
    for date_str in sorted(full_median_curve.keys()):
        med_val = full_median_curve[date_str]
        if med_val > median_peak_swe:
            median_peak_swe = med_val
            median_peak_date = date_str

    if not median_peak_date or median_peak_date <= TODAY.isoformat():
        # If peak hasn't been found in future, something is off — use fallback
        median_peak_date = f"{CURRENT_WY}-04-08"
        # Estimate peak from historical data
        peak_doy = date_to_doy(median_peak_date, CURRENT_WY - 1)
        vals = [historical_curves[yr].get(peak_doy, 0) for yr in historical_curves]
        vals = [v for v in vals if v > 0]
        if vals:
            vals.sort()
            median_peak_swe = vals[len(vals) // 2]
        else:
            median_peak_swe = median_today * 1.5

    print(f"  Median peak: {median_peak_swe} in on {median_peak_date}")

    peak_dt = datetime.date.fromisoformat(median_peak_date)
    days_to_peak = max(1, (peak_dt - TODAY).days)
    swe_needed = round(median_peak_swe - current_swe, 2)
    daily_rate_needed = round(swe_needed / days_to_peak, 4) if days_to_peak > 0 else 0

    # Actual 30-day accumulation rate
    date_30d_ago = (TODAY - datetime.timedelta(days=30)).isoformat()
    swe_30d_ago = None
    for d in sorted_dates:
        if d <= date_30d_ago:
            swe_30d_ago = current_agg[d]["swe"]
    if swe_30d_ago is not None:
        actual_rate = round((current_swe - swe_30d_ago) / 30, 4)
        recent_accum = round(current_swe - swe_30d_ago, 2)
    else:
        actual_rate = 0
        recent_accum = 0

    # Normal remaining accumulation (what typically falls between now and peak)
    median_at_peak = full_median_curve.get(median_peak_date, median_peak_swe)
    normal_remaining = round(median_at_peak - median_today, 2)
    pct_normal_needed = round((swe_needed / normal_remaining) * 100, 1) if normal_remaining > 0 else 999

    # Generate catch-up projection line
    catchup_dates = []
    catchup_swe = []
    if days_to_peak > 0:
        for i in range(days_to_peak + 1):
            d = TODAY + datetime.timedelta(days=i)
            catchup_dates.append(d.isoformat())
            catchup_swe.append(round(current_swe + daily_rate_needed * i, 2))

    # -------------------------------------------------------
    # Step 7: Per-basin stats
    # -------------------------------------------------------
    print("Computing per-basin stats...")
    basin_stats = {}
    for basin_name, triplets in sorted(basin_triplets.items()):
        basin_agg = aggregate_by_date(current_data, station_filter=triplets)
        if basin_agg:
            latest = sorted(basin_agg.keys())[-1]
            b_swe = basin_agg[latest]["swe"]
            b_med = basin_agg[latest]["median"]
            basin_stats[basin_name] = {
                "station_count": len(triplets),
                "current_swe": b_swe,
                "median_swe_today": b_med,
                "pct_of_median": round((b_swe / b_med) * 100, 1) if b_med > 0 else 0,
                "deficit_inches": round(b_med - b_swe, 2),
            }

    # -------------------------------------------------------
    # Step 8: Prepare historical year curves for output
    # Map all historical years to current WY dates for chart overlay
    # -------------------------------------------------------
    hist_output = {}
    for yr, curve in historical_curves.items():
        sorted_doys = sorted(curve.keys())
        # Map DOY to current water year calendar for chart alignment
        hist_output[yr] = {
            "dates": [doy_to_date_str(doy, CURRENT_WY - 1) for doy in sorted_doys],
            "swe": [curve[doy] for doy in sorted_doys],
        }

    # -------------------------------------------------------
    # Step 9: Assemble and write output
    # -------------------------------------------------------
    print("\nWriting data.js...")

    output = {
        "generated_at": datetime.datetime.now().isoformat(timespec="seconds"),
        "water_year": CURRENT_WY,
        "stations": {
            "count": len(stations),
            "list": [
                {
                    "triplet": s["triplet"],
                    "name": s["name"],
                    "elevation": s["elevation"],
                    "lat": s["lat"],
                    "lon": s["lon"],
                    "basin": s["basin"],
                }
                for s in stations
            ],
        },
        "current_season": {
            "dates": sorted_dates,
            "swe": [current_agg[d]["swe"] for d in sorted_dates],
            "median": [current_agg[d]["median"] for d in sorted_dates],
        },
        "median_full_season": {
            "dates": sorted(full_median_curve.keys()),
            "median": [full_median_curve[d] for d in sorted(full_median_curve.keys())],
        },
        "historical_envelope": {
            "dates": date_labels,
            "min": envelope["min"],
            "p10": envelope["p10"],
            "p25": envelope["p25"],
            "p75": envelope["p75"],
            "p90": envelope["p90"],
            "max": envelope["max"],
        },
        "historical_years": hist_output,
        "summary": {
            "today": TODAY.isoformat(),
            "current_swe": current_swe,
            "median_swe_today": median_today,
            "deficit_inches": deficit,
            "pct_of_median": pct_of_median,
            "median_peak_swe": round(median_peak_swe, 2),
            "median_peak_date": median_peak_date,
            "days_to_peak": days_to_peak,
            "swe_needed_for_peak": swe_needed,
            "daily_rate_needed": daily_rate_needed,
            "actual_daily_rate_30d": actual_rate,
            "pct_of_normal_remaining_needed": pct_normal_needed,
            "recent_accumulation_30d": recent_accum,
        },
        "catchup_projection": {
            "dates": catchup_dates,
            "swe": catchup_swe,
        },
        "basins": basin_stats,
        "forecast_points": [
            {"name": "Steamboat / Yampa", "lat": 40.48, "lon": -106.83},
            {"name": "Vail / Eagle", "lat": 39.64, "lon": -106.37},
            {"name": "Aspen / Gunnison", "lat": 39.19, "lon": -106.82},
            {"name": "Wolf Creek / San Juan", "lat": 37.47, "lon": -106.80},
        ],
    }

    # Write as JavaScript variable for easy loading in HTML
    json_str = json.dumps(output, indent=2)
    with open("data.js", "w") as f:
        f.write(f"const DATA = {json_str};\n")

    # Also write raw JSON for reference
    with open("data.json", "w") as f:
        f.write(json_str)

    print(f"\nDone! Files written:")
    print(f"  data.js  ({len(json_str) // 1024} KB)")
    print(f"  data.json ({len(json_str) // 1024} KB)")
    print(f"\nSummary:")
    print(f"  Stations: {len(stations)}")
    print(f"  Current SWE: {current_swe} in")
    print(f"  Median for today: {median_today} in")
    print(f"  % of median: {pct_of_median}%")
    print(f"  Deficit: {deficit} in")
    print(f"  Peak target: {median_peak_swe} in by {median_peak_date}")
    print(f"  Daily rate needed: {daily_rate_needed} in/day")
    print(f"  Actual 30-day rate: {actual_rate} in/day")


if __name__ == "__main__":
    main()
