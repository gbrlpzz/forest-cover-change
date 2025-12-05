# Vegetation Cover Change Detection via NDVI Trend Analysis on GEE

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17831391.svg)](https://doi.org/10.5281/zenodo.17831391)

A Google Earth Engine script for detecting multi-decadal vegetation change using Landsat Collection 2 time series.

## Quick Start

1. Open [Google Earth Engine Code Editor](https://code.earthengine.google.com/)
2. Copy `ndvi_threshold_change.js` contents
3. Click **Run**
4. Click anywhere on the map to inspect points

## Features

- **8 Change Classes**: Loss, Degradation, Emerging, Maturation, Densification, Establishment, Sparse Accumulation, and Transitional Accumulation
- **Dynamic Configuration**: Easily adjustable analysis period (Year Start/End) and thresholds
- **Hemispheric/Seasonal Adaptation**: Configurable START/END months for Northern/Southern hemisphere analysis
- **Sensitivity Analysis**: Built-in parameter to test threshold stability
- **Multi-scale Trend Analysis**: Compares long-term trend with recent short-term trend
- **Trend Acceleration**: Identifies accelerating vs decelerating growth in the Inspector
- **Statistical Significance**: Trends filtered by Mann-Kendall test (p < 0.05) to reject noise.
- **Collection 2 Natively**: Uses USGS Collection 2 Level-2 Surface Reflectance directly without legacy harmonization.
- **Dynamic Legend**: Updates based on active layer and time configuration
- **Point Inspector**: NDVI, trend, classification, and projection
- **Epoch Tracking**: When areas first reached dense canopy (dynamically generated epochs)
- **Trajectory Projection**: Estimated year to reach dense canopy

## Layers

| Layer | Default | Description |
|-------|---------|-------------|
| Vegetation Change | ✓ On | 8-class change detection |
| Canopy Establishment Epoch | Off | When area first reached dense canopy |
| Years to Dense Canopy | Off | Projected years for gaining areas |

## Classification

| NDVI Class | Range | Trend Class | Slope |
|------------|-------|-------------|-------|
| Dense | ≥ 0.6 | Gaining | > +0.005/yr |
| Transitional | 0.4-0.6 | Stable | ±0.005/yr |
| Sparse | 0.2-0.4 | Losing | < -0.005/yr |
| Bare | < 0.2 | | |

## Limitations

- **Thresholds are approximate**: Optimal values vary by region and ecosystem
- **Linearity Assumption**: The "Years to Dense Canopy" projection is a theoretical signal, not an ecological prediction.
- **Validation Status (Visual Only)**: This tool is experimental. Accuracy has been assessed visually but not quantitatively.
- **Sensor Homogeneity**: Minor spectral differences (TM vs OLI) are uncorrected but deemed acceptable for Collection 2.
- **30m resolution**: May not capture fine-scale patterns

See [docs/methodology.pdf](docs/methodology.pdf) for documentation, limitations, and references.

## Data Sources

- Landsat 5/7/8/9 Surface Reflectance (Collection 2, Tier 1)
- 1985–2025 analysis period

## Citation

If you use this software in your research, please cite:

**APA:**
```
Pizzi, G. (2025). Vegetation Cover Change Detection via NDVI Trend Analysis on GEE (Version 2.0.0) [Computer software]. https://doi.org/10.5281/zenodo.17831391
https://github.com/gbrlpzz/ndvi-vegetation-cover-change
```

**BibTeX:**
```bibtex
@software{pizzi2025vegetation,
  author       = {Pizzi, Gabriele},
  title        = {Vegetation Cover Change Detection via NDVI Trend Analysis on GEE},
  year         = 2025,
  version      = {2.0.0},
  doi          = {10.5281/zenodo.17831391},
  url          = {https://github.com/gbrlpzz/ndvi-vegetation-cover-change}
}
```

For the methodological documentation, see [docs/methodology.pdf](docs/methodology.pdf).

## License

Apache License 2.0

## Contact

Gabriele Pizzi | info@gabrielepizzi.com
