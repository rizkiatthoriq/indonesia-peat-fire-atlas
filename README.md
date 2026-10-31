# Indonesia peatland fire activity dashboard

A Google Earth Engine Code Editor dashboard for exploring where high-confidence MODIS active-fire detection pixels occur on mapped Indonesian peatland. It compares individual months with a cumulative Gaussian-smoothed density surface from October 2025 through October 2026, and places rainfall, burned area, land cover, and nearby residential population in context.

This is a **descriptive portfolio analysis**. It does not predict fire, trace the path of an individual fire, measure burn severity, or estimate smoke exposure.

## Run it

1. Sign in to the [Google Earth Engine Code Editor](https://code.earthengine.google.com/).
2. Create a new script and paste the complete contents of [`indonesia_peat_fire_dashboard.js`](indonesia_peat_fire_dashboard.js).
3. Press **Run**. The default selection is Indonesia; no drawn geometry or private asset is needed.
4. Use the province and peat-class menus to narrow the analysis. Switch between a single month and cumulative activity, or compare 2, 5, and 10 km smoothing distances.
5. Use the map-key checkboxes to view the context layers. Pan or zoom, then press **Update GIF from current view** to render a 13-frame cumulative animation for the visible extent.

Earth Engine access and the public/community catalog assets referenced in the script are required. National calculations and GIF rendering can take several minutes or be limited by the Earth Engine project's compute quota. No raster data are stored in this repository.

## What is shown

| Component | Meaning |
|---|---|
| Fire detections | Daily NASA FIRMS MODIS raster pixels with `confidence > 95`, intersected with mapped peat. A detection pixel can be observed more than once. |
| Monthly and cumulative density | A Gaussian kernel applied to the approximately 1 km detection grid. The default bandwidth is 5 km, with a 15 km kernel radius. Colour breaks stay fixed across months. |
| Peat classes | Global Peatland Map 2.0 pixels classified as peat dominated or peat in a soil mosaic. |
| Burned area | A separate MODIS MCD64A1 burned-pixel estimate, intersected with mapped peat. Its latest available month can lag FIRMS. |
| Rainfall | CHIRPS daily rainfall, summed within 14-day bins and averaged over mapped peat for the chart. The optional map displays rainfall accumulated over available dates. |
| Land cover | Esri/Impact Observatory 2024 annual classes used as pre-period context for areas inside the high density class and mapped burned area. |
| Population | WorldPop 2025 residential counts. Source cells are **summed**, not averaged, onto the approximately 1 km fire grid for the map and distance-band summaries. |
| Nearby residents | Non-overlapping 0–2, 2–5, and 5–10 km ranges from detected peat-fire pixels. They describe proximity, not confirmed impact. |

The code's end date is **30 October 2026 UTC**. That is the analysis window, not a guarantee that every dataset has observations through that day. Check the **latest source dates** displayed in the dashboard before interpreting a result. The October 2026 animation frame covers available observations through the configured end date.

The GIF uses an Earth Engine terrain image and FAO province outlines. The interactive Google Maps basemap cannot be exported as part of the Earth Engine animation.

## Read the results carefully

- `confidence > 95` is a MODIS detection-confidence score, not a statistical 95% confidence interval. NASA describes the Earth Engine FIRMS raster as near-real-time data that are not science quality. [FIRMS catalog](https://developers.google.com/earth-engine/datasets/catalog/FIRMS)
- The density map is a **gridded approximation** to kernel density estimation. It shows repeated detection-pixel concentration; it does not show the number of independent fires or their direction of travel. The 5 km bandwidth and the highest-class threshold of 0.5 approximate observations/km² are exploratory choices, with alternatives exposed in the dashboard. [Diggle (1985)](https://doi.org/10.2307/2347366)
- “Peat in a soil mosaic” gives the area of classified **pixels**, not the area of peat within each mixed pixel. The source map is approximately 1 km. [Global Peatland Map 2.0](https://gee-community-catalog.org/projects/peatland/)
- Density-class overlap in hectares is not burned area. MCD64A1 supplies an independent burned-area layer, with its own date and detection limitations. [MCD64A1 catalog](https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MCD64A1)
- Residents within 10 km are people whose **modeled home locations** fall near a detected pixel. They are not a count of people injured, evacuated, or exposed to smoke. [WorldPop dataset](https://gee-community-catalog.org/projects/worldpop/)
- CHIRPS rainfall is context, not peat-water-table depth or proof that rainfall caused a change in fire detections. [CHIRPS catalog](https://developers.google.com/earth-engine/datasets/catalog/UCSB-CHG_CHIRPS_DAILY)

See [METHODS.md](METHODS.md) for processing steps, units, and attribution.

## Data use and attribution

The repository contains code, not copies of the source datasets. Global Peatland Map 2.0 is supplied under **CC BY-NC-SA 4.0**; this project was prepared for non-commercial portfolio use. Other datasets have their own attribution terms, listed in [METHODS.md](METHODS.md). Check those terms before reusing maps or figures. No software license has been selected for this repository.
