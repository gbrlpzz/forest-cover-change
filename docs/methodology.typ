#set document(
  title: "Vegetation Cover Change Detection via NDVI Trend Analysis: Methodology",
  author: "Gabriele Pizzi",
)
#set page(
  paper: "a4",
  margin: (x: 2.5cm, y: 3cm),
  header: context {
    if counter(page).get().first() > 1 [
      #set text(style: "italic", size: 9pt, fill: gray)
      Vegetation Cover Change Detection via NDVI Trend Analysis | Pizzi 2025
      #h(1fr)
      Methodology
      #line(length: 100%, stroke: 0.5pt + gray)
    ]
  },
  footer: context [
    #line(length: 100%, stroke: 0.5pt + gray)
    #set text(size: 9pt, fill: gray)
    #align(center)[Page #counter(page).display()]
  ],
)
#set text(font: "New Computer Modern", size: 11pt)
#set heading(numbering: "1.1")
#show heading: it => {
  set text(weight: "bold")
  v(0.5em)
  it
  v(0.3em)
}
#set par(justify: true, leading: 0.65em)
#show figure.where(kind: table): set figure.caption(position: top)
#show figure.caption: set text(size: 10pt)
#show table: set text(size: 10pt)


#align(center + horizon)[
  #text(size: 22pt, weight: "bold")[Vegetation Cover Change Detection via NDVI Trend Analysis]
  #v(0.5cm)
  #text(size: 14pt, style: "italic")[Methodology]
  #v(2cm)
  #text(size: 12pt)[*Gabriele Pizzi*]
  #v(0.5cm)
  #text(size: 11pt)[December 2025]
  
  #v(1fr)
  #text(size: 10pt, fill: gray)[
    © 2025 Gabriele Pizzi. All rights reserved.\
    Version 2.0.0
  ]
]

#pagebreak()

= Abstract

This document details the scientific methodology for detecting multi-decadal vegetation cover change using inter-calibrated Landsat time series. The classification system employs NDVI (Normalized Difference Vegetation Index) thresholds derived from USGS reference values and linear trend analysis parameters adapted from peer-reviewed literature (e.g., Peng & Gong, 2025). The implementation allows for dynamic temporal configuration and sensitivity analysis to adapt to different biomes.

= Introduction

Vegetation monitoring requires consistent definition of land cover states. This methodology implements a standardized taxonomy based on biophysical thresholds to classify vegetation density and change dynamics over a 40-year period.

== Applicability Scope

This methodology is calibrated for temperate and Mediterranean ecosystems. Users should note that optimal NDVI thresholds vary by latitude and biome (Pettorelli et al., 2005). Additionally, the analysis period is configurable to account for hemispheric seasonality differences (e.g., June-September for Northern Hemisphere, December-March for Southern).

= Data Processing

== Satellite Imagery

Collection 2, Level-2 Surface Reflectance data from USGS:

- Landsat 5 TM (1984–2012)
- Landsat 7 ETM+ (1999–2025)
- Landsat 8 OLI (2013–2025)
- Landsat 9 OLI-2 (2021–2025)

All imagery is processed at native 30-meter spatial resolution.

== Quality Masking

Cloud and cloud shadow contamination is removed using the QA_PIXEL quality assessment band following USGS Collection 2 specifications:

- Bit 3 (Cloud) = 0 (clear sky)
- Bit 4 (Cloud Shadow) = 0 (no shadow)

Only clear observations are retained for NDVI calculation, ensuring temporal composites represent actual vegetation conditions rather than atmospheric artifacts.

== Harmonization

== Spectral Harmonization

Landsat 8/9 OLI spectral bands are natively compatible with Landsat 5/7 TM/ETM+ for vegetation analysis when using USGS Collection 2 Level-2 Surface Reflectance data. *Crawford et al. (2023)* demonstrate that the geometric and radiometric improvements in Collection 2 significantly reduce the need for post-hoc OLS harmonization (e.g., Roy et al., 2016) for general monitoring applications.

*Design Decision:* This methodology relies directly on USGS Collection 2 Level-2 inter-calibration. Minor spectral response differences between sensors (OLI vs TM/ETM+) are acknowledged as a known uncertainty, but the preservation of original radiometric data is prioritized over the introduction of potential transformation artifacts.

*Processing Protocol:*
- *Landsat 5 & 7*: Bands B3 (Red) and B4 (NIR).
- *Landsat 8 & 9*: Bands B4 (Red) and B5 (NIR).

Explicit band mapping ensures correct spectral matching without altering pixel values.

= NDVI Classification Thresholds

The core classification relies on absolute NDVI values to define vegetation states.

== Threshold Verification

The selected thresholds align with U.S. Geological Survey (USGS) standards for Remote Sensing Phenology:

#figure(
  table(
    columns: (auto, auto, auto),
    inset: 8pt,
    align: left,
    table.header([*Class*], [*NDVI Range*], [*Standard Interpretation*]),
    table.hline(),
    [Dense Canopy], [≥ 0.6], [Dense forest, healthy vegetation],
    [Transitional], [0.4 – 0.6], [Open woodland, shrubland],
    [Sparse], [0.2 – 0.4], [Grassland, senescing crops],
    [Bare], [< 0.2], [Soil, rock, snow, water],
    table.hline(),
  ),
  caption: [NDVI Classification Thresholds],
)

=== Dense Canopy (≥ 0.6)
Studies on land surface emissivity (Sobrino et al., 2004) classify pixels with NDVI > 0.5 as "fully vegetated." This methodology applies a conservative threshold of 0.6 to define "Dense Canopy," ensuring that only high-biomass, healthy forest structures are captured, significantly reducing false positives from mixed pixels.

=== Sparse Vegetation (0.2 – 0.4)
Sobrino et al. (2004) characterize the range 0.2 ≤ NDVI ≤ 0.5 as "mixed pixels" containing a combination of soil and vegetation. The "Sparse" class (0.2–0.4) strictly targets this heterogeneous interface, while the "Transitional" class (0.4–0.6) captures the upper bound of this mixed zone.

== Sensitivity Analysis

The system includes a `SENSITIVITY_ADJUSTMENT` parameter (default 0.0) that applies a global offset to all NDVI thresholds. This allows researchers to test the stability of classification results against threshold variations (e.g., ±0.05), providing a mechanism to assess the robustness of the detected changes.

#block(
  fill: rgb("#FFF3CD"),
  inset: 10pt,
  radius: 4pt,
  [
    *Design Decision:* The subdivision at 0.4 NDVI is a methodological choice to separate sparse from transitional vegetation classes. This specific threshold is not independently validated in peer-reviewed literature and may require regional calibration for optimal performance in non-temperate ecosystems.
  ]
)

= Trend Analysis

== Linear Trend Computation

Trends are calculated using Ordinary Least Squares (OLS) regression on annual summer median composites. To ensure statistical rigor, the Mann-Kendall trend test (Kendall's Tau) is applied to corresponding pixels.

$ S = sum_(k=1)^(n-1) sum_(j=k+1)^n "sgn"(x_j - x_k) $

Only trends with a statistical significance of $p < 0.05$ (95% confidence level) are retained. Pixels failing this test are classified as "Stable" regardless of their linear slope magnitude.

== Trend Significance Thresholds

To separate natural variability from significant change, a slope threshold of ±0.005 NDVI/year is applied to statistically significant pixels.

#figure(
  table(
    columns: (auto, auto),
    inset: 8pt,
    align: left,
    table.header([*Trend Class*], [*Slope Threshold*]),
    table.hline(),
    [Gaining], [> +0.005/yr AND $p < 0.05$],
    [Stable], [±0.005/yr OR $p >= 0.05$],
    [Losing], [< -0.005/yr AND $p < 0.05$],
    table.hline(),
  ),
  caption: [Trend Significance Thresholds],
)

=== Scientific Basis
This threshold is derived from Peng & Gong (2025), whose analysis of spatiotemporal NDVI changes classified slopes between 0.005 and 0.016 as "moderate improvement" and defined the stable range as -0.007 to 0.005. This provides a peer-reviewed basis for the cutoff.

=== Recent Trend Analysis

To detect acceleration or deceleration in recent vegetation change, a secondary 10-year trend (2015–2025) is computed and compared to the 40-year baseline trend. This trend acceleration indicator identifies whether change is intensifying or moderating:

#figure(
  table(
    columns: (auto, auto),
    inset: 8pt,
    align: left,
    table.header([*Change Rate Class*], [*Condition*]),
    table.hline(),
    [Accelerating], [Recent slope > Long-term slope + 0.002],
    [Consistent], [|Recent slope - Long-term slope| < 0.002],
    [Decelerating], [Recent slope < Long-term slope - 0.002],
    table.hline(),
  ),
  caption: [Recent Trend Classification Criteria],
)

The 0.002 NDVI/year threshold was selected to distinguish meaningful acceleration from noise while remaining sensitive to ecological change dynamics. This dual-timeframe approach helps identify recent shifts in land management or climate-driven vegetation responses.

= Classification Taxonomy

The system intersects absolute state (NDVI) with directional trend (Slope) to produce 8 mutually exclusive classes. To ensure robustness against statistical noise, the classification prioritizes *State Change* (difference between start/end median NDVI) over linear trends. Linear trends are used as a secondary confirmation for subtle intra-class changes (e.g., Densification, Accumulation).

== State-Driven Classes (Transition Logic)

#figure(
  table(
    columns: (auto, auto, auto),
    inset: 8pt,
    align: left,
    table.header([*Class*], [*Definition*], [*Transition Logic*]),
    table.hline(),
    [Canopy Loss], [Dense → Sparse/Bare], [State collapse],
    [Degradation], [Dense → Trans], [Biomass decline],
    [Emerging Biomass], [Sparse → Trans], [Early recovery phase],
    [Maturation], [Trans → Dense], [Full canopy closure],
    [Establishment], [Sparse/Bare → Dense], [Rapid afforestation],
    [Densification], [Dense → Dense (+Gain)], [Biomass increase (Trend-driven)],
    [Transitional Accumulation], [Trans → Trans (+Gain)], [Intra-class growth (approaching Dense)],
    [Sparse Accumulation], [Sparse → Sparse (+Gain)], [Intra-class growth (approaching Trans)],
    table.hline(),
  ),
  caption: [Simplified Classification Matrix],
)

== Canopy Establishment Epochs

For areas classified as "Canopy Establishment" (Sparse/Bare → Dense), the specific time period when dense canopy was first achieved is tracked using 5-year epochs generated dynamically based on the analysis period.

#figure(
  table(
    columns: (auto, auto, auto),
    inset: 8pt,
    align: center,
    table.header([*Epoch Label*], [*Time Period*], [*Interpretation*]),
    table.hline(),
    [1990], [1990–1994], [Early establishment],
    [1995], [1995–1999], [],
    [2000], [2000–2004], [Millennium era],
    [2005], [2005–2009], [],
    [2010], [2010–2014], [Recent decade],
    [2015], [2015–2019], [],
    [2020], [2020–2025], [Latest period],
    table.hline(),
  ),
  caption: [Canopy Establishment Epoch Definitions],
)

The baseline period (1985–1989) is excluded from epoch tracking as it serves as the initial reference state. The 5-year interval balances temporal precision with data availability, ensuring sufficient cloud-free observations for robust NDVI composites within each epoch.

// Edge classes removed in version 2.0.0 due to stable trend masking.

= Limitations & Caveats

Users must acknowledge the following limitations when interpreting results:

1.  *Threshold Universality*: While the >0.5 threshold for full vegetation is standard (Sobrino et al., 2004), the use of 0.6 is conservative. Optimal values vary by region, and boreal or dryland forests may require lower thresholds.
2.  *Validation Status (Visual Only)*: This methodology has undergone *preliminary visual validation* by the author in select test sites. It has *NOT* been rigorously validated with a quantitative accuracy assessment (confusion matrix). Users should treat the output as experimental indices rather than ground-truth maps until further validation is published.
3.  *Linearity Assumption*: The "Years to Dense Canopy" projection is a *theoretical linear model* (`(Threshold - Current) / RecentSlope`). It projects the *current 10-year trend* forward. Ecological recovery is typically sigmoid/asymptotic. This projection likely underestimates recovery time for young stands (slow start) and typically fails to account for carrying capacity saturation in mature stands.
4.  *Sensor Homogeneity*: Despite Collection 2 inter-calibration (Crawford et al., 2023), minor spectral differences between Landsat generations (TM/ETM+ vs OLI) may persist. It is assumed these are negligible for broad degradation classes, but they may influence subtle trend detection across the 2012/2013 sensor transition.
5.  *Projection Limits*: The "Years to Dense Canopy" layer is capped at 50 years (`clamp(0, 50)`) for visualization purposes. This artificial horizon should be considered when interpreting long-term recovery projections.

= Code Availability

The complete source code, including the Google Earth Engine script and documentation, is available at:

#link("https://github.com/gbrlpzz/ndvi-vegetation-cover-change")

This implementation was developed using the Google Earth Engine JavaScript API via the Code Editor interface (tested as of December 2025). While GEE maintains backward compatibility, users should be aware that API updates may occasionally require minor code adjustments for future compatibility.

= References

#set par(hanging-indent: 2em)

Crawford, C. J., et al. (2023). The 50-year Landsat collection 2 archive. _Science of Remote Sensing_, 8, 100103. #link("https://doi.org/10.1016/j.srs.2023.100103")
 
Peng, Y., & Gong, H. (2025). Analysis of Spatiotemporal Changes in NDVI-Derived Vegetation Index and Its Influencing Factors in Kunming City (2000 to 2020). _Forests_, 16(12), 1781. #link("https://doi.org/10.3390/f16121781")
 
Pettorelli, N., Vik, J. O., Mysterud, A., Gaillard, J. M., Tucker, C. J., & Stenseth, N. C. (2005). Using the satellite-derived NDVI to assess ecological responses to environmental change. _Trends in Ecology & Evolution_, 20(9), 503–510. #link("https://doi.org/10.1016/j.tree.2005.05.011")

Sobrino, J. A., Jiménez-Muñoz, J. C., & Paolini, L. (2004). Land surface temperature retrieval from LANDSAT TM 5. _Remote Sensing of Environment_, 90(4), 434–440. #link("https://doi.org/10.1016/j.rse.2004.02.003")

#v(0.5cm)
#align(center)[
  #text(size: 9pt, fill: gray)[
    Document version 2.0.0 | December 2025 | © 2025 Gabriele Pizzi
  ]
]
