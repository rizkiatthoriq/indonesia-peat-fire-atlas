# Methods, units, and source credits

## Analysis settings

| Setting | Value in the script | Interpretation |
|---|---|---|
| Study area | Indonesia, FAO GAUL 2025 level 1; province menu limited to units overlapping mapped peat | National maps include all Indonesian provinces; fire analysis uses mapped peat pixels within them. |
| Time window | 1 October 2025 to 30 October 2026 UTC, inclusive | Earth Engine uses `2026-10-31` as the exclusive filter end. Individual source collections may end earlier. |
| Active-fire filter | FIRMS MODIS `confidence > 95` | A quality-screening score on rasterized active-fire detections, not a confidence interval. |
| Gaussian smoothing | 5 km sigma, 15 km kernel radius by default; 2 and 10 km sigma options | Descriptive smoothing of a roughly 1 km detection grid. The radius is three times sigma. Bandwidth has not been optimized statistically. |
| Density display | Fixed class starts: 0.001, 0.01, 0.05, 0.15, 0.5 approximate detection observations/km² | Fixed breaks support visual month-to-month comparison. The highest class is an exploratory display threshold. |
| Chart interval | 14 days | Bars count detection-pixel observations; the rainfall line is the spatial mean of each bin's CHIRPS total over mapped peat, in mm. |
| Population ranges | 0–2, 2–5, and 5–10 km | Exclusive bands from the nearest detected peat-fire pixel. The 0–10 km card adds the three bands. |

## Processing sequence

1. Select GAUL level 1 units for Indonesia and the peat-dominated and peat-mosaic pixels of Global Peatland Map 2.0.
2. Filter daily Earth Engine FIRMS MODIS rasters to `confidence > 95`. Count qualifying **pixels per daily image**, then sum across the chosen month or all months to date. Repeated satellite observations can count repeatedly.
3. Convolve the detection-count raster with a normalized Gaussian kernel. Divide by pixel area in km² to obtain an approximate observation intensity per km². Restrict the displayed density to mapped peat and apply fixed colour breaks. This is a raster-grid KDE approximation, not an analysis of original vector point coordinates, and the script does not implement an edge correction.
4. For each successive month, recompute from 1 October 2025 through that month's end. The cumulative animation therefore shows changing **concentration**, not the movement of an individual fire front.
5. Separately intersect valid MCD64A1 burned pixels with mapped peat. Report mapped burned hectares and their approximate overlap with the highest density class. The high-density footprint itself is also reported in hectares, but is not called burned area.
6. Use 2024 Esri/Impact Observatory land cover as pre-period context. The national high-density class summary samples at 1 km for computational feasibility; province summaries use 100 m. Treat resulting class hectares as estimates.
7. Sum WorldPop's modeled 2025 residential **people per cell** from its approximately 100 m grid to the FIRMS grid. Measure grid-cell distance to the nearest qualifying peat-fire pixel and sum residents in exclusive 0–2, 2–5, and 5–10 km ranges. Residents can live off peat. These ranges are descriptive proximity choices, not validated exposure thresholds.
8. Sum CHIRPS daily rainfall within each 14-day interval, then average the interval total spatially across mapped peat. Missing or incomplete rainfall bins are left blank. A failed bin is reported in the GEE Console rather than silently treated as zero.

## Data sources and reuse conditions

| Layer | Source and credit | Resolution / key condition |
|---|---|---|
| Peat | [Greifswald Mire Centre (2022), Global Peatland Map 2.0](https://gee-community-catalog.org/projects/peatland/), curated in the GEE Community Catalog | ~1 km; classes 1 and 2; **CC BY-NC-SA 4.0**. |
| Administrative boundaries | [FAO GAUL 2025 level 1](https://developers.google.com/earth-engine/datasets/catalog/FAO_GAUL_2025_level1) | Province-level polygons; **CC BY 4.0**. |
| Active fire | [NASA LANCE/FIRMS MODIS NRT](https://developers.google.com/earth-engine/datasets/catalog/FIRMS), DOI [10.5067/FIRMS/MODIS/MCD14DL.NRT.006](https://doi.org/10.5067/FIRMS/MODIS/MCD14DL.NRT.006) | ~1 km rasterized detections; source requests citation and disclaimer when data are passed to others. |
| Rainfall | [UCSB CHG CHIRPS daily](https://developers.google.com/earth-engine/datasets/catalog/UCSB-CHG_CHIRPS_DAILY); [Funk et al. (2015)](https://doi.org/10.1038/sdata.2015.66) | ~5.6 km; daily mm; public domain. |
| Burned area | [NASA MODIS MCD64A1 Collection 6.1](https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MCD64A1); [Giglio et al. (2018)](https://doi.org/10.1016/j.rse.2018.08.005) | 500 m monthly burned-pixel and QA bands. |
| Land cover | [Impact Observatory for Esri, annual land cover](https://gee-community-catalog.org/projects/S2TSLULC/); [Karra et al. (2021)](https://doi.org/10.1109/IGARSS47720.2021.9553499) | 2024 source at 10 m; **CC BY 4.0**. Credit: “Produced by Impact Observatory for Esri. © 2021 Esri.” |
| Population | [WorldPop Global 2015–2030](https://gee-community-catalog.org/projects/worldpop/), total population model | 2025 constrained residential counts at ~100 m; **CC BY 4.0**; summed to ~1 km in this script. |

The GIF also uses [NOAA ETOPO1](https://developers.google.com/earth-engine/datasets/catalog/NOAA_NGDC_ETOPO1) for an exportable terrain base and `users/gena/packages:text` for month labels. Review the terms of these dependencies before redistributing a rendered animation.

## Validation and limits

The JavaScript file is syntax checked locally. Earlier interactive GEE checks confirmed the nationwide rainfall overlay, peatland province filter, and national and Central Kalimantan population-band calculations. Those checks used the same analytical code before the latest panel-size and boundary-width adjustments. The final visual map state has **not** been fully revalidated after the Earth Engine project displayed a non-commercial compute-quota restriction. Source coverage and output values should be checked in a fresh run before publishing numerical results or screenshots.

The dashboard reads source dates at runtime. It should not be described as complete through 30 October 2026 unless the relevant source actually contains that date. The near-real-time FIRMS collection is not science quality; the peat and land-cover classifications and modeled population each add uncertainty. The analysis does not estimate ignition cause, spread direction, smoke concentration, exposure, or ecological loss.
