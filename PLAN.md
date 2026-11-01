# Research plan

## Question
How has snow cover changed since 2000 across the Upper Indus and the wider Hindu Kush-Karakoram-Himalaya, and which climate signals explain the differences between basins, ranges and elevations?

## Regions, nested and reported separately
Upper Indus and its sub-basins; Third Pole domain; Hindu Kush, Karakoram, Himalaya (west, central, east) and Pamir. Boundaries are fixed once, versioned, and defined by hydrology and topography, not national borders. Each is also split into elevation bands.

## Data (checked in Earth Engine, 2026-11-01)
- Snow: MODIS/061/MOD10A1 and MYD10A1 (daily), MOD10A2 (8-day). MOD10A1F, the cloud-gap-filled product, is not in the catalogue: test gap filling ourselves, or process NSIDC files locally.
- Climate: ERA5-Land monthly, IMERG V07 monthly, CHIRPS (up to 50N, limited use at altitude), OISST V2.1.
- Terrain and ice: Copernicus DEM, HydroBASINS, GLIMS or RGI.
- Independent checks: Landsat and Sentinel-2 snow maps; station records where obtainable.

## Metrics, per hydrological year
Snow cover area: mean, median, min, max; snow-cover duration; onset and melt-out dates; snowline elevation; each with the fraction of valid observations.

## Method
1. Reproduce `original/UIBNDSI.xlsx` with MOD10A2 on the original domain to establish continuity.
2. Handle cloud explicitly, and compare at least two treatments; report the trend's sensitivity to them.
3. Validate against Landsat and Sentinel-2 by elevation band and season.
4. Trends: Mann-Kendall with Sen slope, allowing for autocorrelation, with field significance control across many regions.
5. Attribution: temperature, precipitation and snowfall; ENSO, IOD, western disturbance and monsoon indices with lags, by season. Report associations with uncertainty; do not claim causation from correlation.

## Where it runs
Earth Engine: compositing, zonal statistics, exports to CSV. Local Python: trends, attribution, figures, tables, tests. One script per figure.

## Reuse
`thirdpole3-5` and `monsoon8` in the Earth Engine repository, consolidated into shared functions. Each is tested before use.

## Writing
New text only. Earlier reports and slides are for reference, as they contain material that has not been checked for sources.
