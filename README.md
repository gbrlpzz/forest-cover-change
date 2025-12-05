# Vegetation Cover Change Detection

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

A Google Earth Engine script for detecting 40-year vegetation change using harmonized Landsat time series.

## Quick Start

1. Open [Google Earth Engine Code Editor](https://code.earthengine.google.com/)
2. Copy `ndvi_threshold_change.js` contents
3. Click **Run**
4. Click anywhere on the map to inspect points

## Features

- **9 Change Classes**: Canopy loss, thinning, emerging biomass, thickening, densification, establishment, and 3 edge dynamics
- **Multi-scale Trend Analysis**: Compares 40-year long-term trend with 10-year recent trend (Momentum)
- **Momentum Indicator**: Identifies accelerating vs decelerating growth in the Inspector
- **Spectral Harmonization**: Landsat 8/9 data harmonized to Landsat 7 standards (Roy et al., 2016)
- **Dynamic Legend**: Updates based on active layer
- **Point Inspector**: NDVI, trend, classification, and projection
- **Epoch Tracking**: When areas first reached dense canopy
- **Trajectory Projection**: Estimated year to reach dense canopy

## Layers

| Layer | Default | Description |
|-------|---------|-------------|
| Vegetation Change | ✓ On | 9-class change detection |
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
- **Linear trend assumption**: Vegetation change is often nonlinear
- **No ground validation**: Results should be validated locally
- **30m resolution**: May not capture fine-scale patterns

See [docs/methodology.pdf](docs/methodology.pdf) for documentation, limitations, and references.

## Data Sources

- Landsat 5/7/8/9 Surface Reflectance (Collection 2, Tier 1)
- 1985–2025 analysis period

## Citation

If you use this software in your research, please cite:

**APA:**
```
Pizzi, G. (2025). Vegetation Cover Change Detection (Version 2.0.0) [Computer software]. 
https://doi.org/10.5281/zenodo.XXXXXXX
```

**BibTeX:**
```bibtex
@software{pizzi2025vegetation,
  author       = {Pizzi, Gabriele},
  title        = {Vegetation Cover Change Detection},
  year         = 2025,
  version      = {2.0.0},
  doi          = {10.5281/zenodo.XXXXXXX},
  url          = {https://github.com/gbrlpzz/ndvi-vegetation-cover-change}
}
```

For the methodological documentation, see [docs/methodology.pdf](docs/methodology.pdf).

## License

Apache License 2.0

## Contact

Gabriele Pizzi | info@gabrielepizzi.com
