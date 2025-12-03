# Forest Change Detection (1985-2023) using Landsat

## Overview

This Google Earth Engine (GEE) script performs a retrospective analysis of forest change (loss and gain) between 1985 and 2023 for a specified Region of Interest (ROI). It utilizes harmonized surface reflectance data from Landsat 5, 7, 8, and 9 missions to generate two primary outputs:

Deforestation (Loss) Mask: Areas classified as forest during the baseline period (1985–1987) that transitioned to open land by the end period (2021–2023).

Reforestation (Gain) Map: Areas classified as open land during the baseline period that transitioned to forest by the end period. For these gain areas, the map provides the approximate year of detection, the year the pixel first crossed the forest threshold.

The workflow for regions where abandonment and forest transition occurred before the Sentinel-2 era, requiring multi-decadal Landsat continuity.

The results are displayed on the map with color coding and configured for export as GeoTIFFs to Google Drive.

## Methodology

The analysis follows a standard land cover change-detection workflow, using the Normalized Difference Vegetation Index (NDVI) as the primary proxy for forest cover.

1. Data Harmonization

All available Landsat 5, 7, 8, and 9 Level-2, Collection 2, Tier 1 Surface Reflectance data are merged. Each image is processed to mask clouds and cloud shadows using the QA_PIXEL band. Scale factors and offsets are applied to convert the data to standard surface reflectance values, and bands are renamed to NIR and Red for consistency.

2. Forest Classification

A fixed NDVI threshold of 0.45 (forestThreshold) is used to classify pixels:

NDVI > 0.45: Forest

NDVI ≤ 0.45: Open Land / Non-Forest

3. Change Calculation

Change is determined by comparing the median summer NDVI (June, July, August, September) across two distinct periods:

Baseline (Start State): 1985–1987
Final (End State): 2021–2023

Deforestation (Loss): Pixels that were Forest in the Start State and Open in the End State.
Reforestation (Gain): Pixels that were Open in the Start State and Forest in the End State.

4. Year of Recovery

For pixels identified as Reforestation, a pixel-by-pixel time series analysis is performed between startYear (1985) and endYear (2023). The script finds the earliest year in this range where the pixel's median summer NDVI first crosses the 0.45 forest threshold. This year is recorded in the finalRecoveryMap.

## Example Output

<img width="1289" height="906" alt="Screenshot 2025-12-03 at 18 18 59" src="https://github.com/user-attachments/assets/e2533ce5-cfd6-448f-8f5d-583edd25e035" />

## Intended Applications

This workflow is designed for analyzing long-term land use transitions where forest cover gain or loss reflects broader socio-ecological processes. Typical use-cases include:

* **Agricultural land abandonment** and subsequent **forest encroachment**
* **Forest recovery** following disturbances, logging, or historical land-use changes
* **Landscape transition monitoring** in rural and peri-urban territories
* **Temporal mapping** of legacy land-use decisions with ecological implications

While Sentinel-2 offers higher spatial resolution (10 m), this system uses **harmonized Landsat data (30 m) since 1985** to capture **decades-long trajectories** of abandonment and forest succession. This provides the temporal depth required in regions—such as much of Italy—where major land abandonment and recovery events occurred well before the Copernicus era.

## Key Variables

| Variable | Description | Default Value |
| startYear | Start year for annual change analysis. | 1985 |
| endYear | End year for annual change analysis. | 2023 |
| forestThreshold | NDVI value used to define forest cover. | 0.45 |
| roi | The region of interest for analysis. | Current Map Bounds |

## Data Sources

Landsat 5 TM: LANDSAT/LT05/C02/T1_L2

Landsat 7 ETM+: LANDSAT/LE07/C02/T1_L2

Landsat 8 OLI/TIRS: LANDSAT/LC08/C02/T1_L2

Landsat 9 OLI-2/TIRS-2: LANDSAT/LC09/C02/T1_L2



## Contact

Gabriele Pizzi | info@gabrielepizzi.com | gabrielepizzi.com
