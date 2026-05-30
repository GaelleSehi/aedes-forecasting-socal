# Boundary Shapefiles

This directory contains two shapefiles required to run the spatial analyses.
They are not included in the repository due to file size constraints.

---

## 1. Study-area boundary (`inst/boundaries/district/`)

The district boundary polygon for the San Gabriel Valley Mosquito and Vector
Control District study area.

**Source:** Provided by the San Gabriel Valley Mosquito and Vector Control
District. Available upon request from the corresponding author.

**Expected format:** A single polygon shapefile. The layer name used in the
scripts is `"Untitled_layer-polygon"` — update `boundary_path` in
`R/00_config.R` if your layer name differs.

**CRS:** Any standard geographic CRS is accepted; scripts reproject to
EPSG:4326 (WGS84) and EPSG:32611 (UTM Zone 11N) as needed.

---

## 2. California counties (`inst/boundaries/Cali_3310/`)

California county boundaries, used as a background layer in error maps.

**Source:** California Department of Technology / State of California Open Data.
Download from: https://data.ca.gov/dataset/ca-geographic-boundaries

- Direct download (counties): `ca-county-boundaries.zip`
- Or via the `tigris` R package:

```r
library(tigris)
ca <- counties(state = "CA", cb = TRUE, year = 2020)
sf::st_write(ca, "inst/boundaries/Cali_3310/CA_counties.shp")
```

**CRS:** Scripts reproject to match the study-area boundary CRS (EPSG:32611).

---

## Notes

- After downloading, place the extracted `.shp`, `.shx`, `.dbf`, and `.prj`
  files into the respective subdirectory.
- Paths are defined in `R/00_config.R` as `boundary_path` and `counties_path`.
- Both directories are listed in `.gitignore`.
