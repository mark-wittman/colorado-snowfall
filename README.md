# Colorado Snowpack Dashboard

Shows how Colorado's snowpack compares to historical norms and what's needed to catch up.

## Quick Start

1. Install Python dependency (one time):
   ```
   cd ~/Projects/colorado-snowfall
   source .venv/bin/activate
   ```

2. Fetch the latest data:
   ```
   python fetch_data.py
   ```

3. Open the dashboard:
   ```
   open index.html
   ```

## How to Update

Run `python fetch_data.py` again whenever you want fresh SNOTEL data. The 7-day
forecast updates automatically each time you open the dashboard (fetched live).

## What's on the Dashboard

- **Hero cards**: Current SWE, deficit, % of median, days to peak
- **Season curve**: This year vs 20-year median, with 10th-90th percentile historical range
- **Catch-up analysis**: Pace needed vs actual pace, rate multiple, % of normal remaining
- **7-day forecast**: Expected mountain snowfall from 4 CO locations (Open-Meteo)
- **Basin breakdown**: Per-basin SWE stats for 6 Colorado river basins
- **Basin conditions**: Per-basin 24h average temperature, wind, and near-surface soil moisture (Open-Meteo)

## Data Sources

- **SNOTEL**: NRCS AWDB REST API — 119 Colorado stations, no auth needed
- **Forecast**: Open-Meteo API — free, no auth needed
