#set document(
  title: "Vegetation Cover Change Detection: Methodology",
  author: "Gabriele Pizzi",
)
#set page(
  paper: "a4",
  margin: (x: 2.5cm, y: 3cm),
  header: context {
    if counter(page).get().first() > 1 [
      #set text(style: "italic", size: 9pt, fill: gray)
      Vegetation Cover Change Detection | Pizzi 2025
      #h(1fr)
      Methodology & Validation
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
  #text(size: 22pt, weight: "bold")[Vegetation Cover Change Detection]
  #v(0.5cm)
  #text(size: 14pt, style: "italic")[Methodology and Scientific Validation]
  #v(2cm)
  #text(size: 12pt)[*Gabriele Pizzi*]
  #v(0.5cm)
  #text(size: 11pt)[December 2025]
  
  #v(1fr)
  #text(size: 10pt, fill: gray)[
    © 2025 Gabriele Pizzi. All rights reserved.\
    Version 2.2
  ]
]

#pagebreak()

= Abstract

This document details the scientific methodology for detecting multi-decadal vegetation cover change (1985–2025) using harmonized Landsat time series. The classification system employs NDVI (Normalized Difference Vegetation Index) thresholds derived from USGS standards and linear trend analysis parameters adapted from peer-reviewed literature (e.g., Peng & Gong, 2025).

= Introduction

Vegetation monitoring requires consistent definition of land cover states. This methodology implements a standardized taxonomy based on biophysical thresholds to classify vegetation density and change dynamics over a 40-year period.

== Applicability Scope

This methodology is calibrated for temperate and Mediterranean ecosystems. Users should note that optimal NDVI thresholds vary by latitude and biome (Pettorelli et al., 2005).

= Data Processing

== Satellite Imagery

Collection 2, Level-2 Surface Reflectance data from USGS:

- Landsat 5 TM (1984–2012)
- Landsat 7 ETM+ (1999–2025)
- Landsat 8 OLI (2013–2025)
- Landsat 9 OLI-2 (2021–2025)

== Harmonization

== Spectral Harmonization

Landsat 8/9 OLI spectral response differs from Landsat 5/7 TM/ETM+. Following Roy et al. (2016), band values are harmonized to ensure temporal consistency using Ordinary Least Squares (OLS) regression coefficients.

#figure(
  table(
    columns: (auto, auto, auto),
    inset: 8pt,
    align: center,
    stroke: none,
    table.header([*Band*], [*Slope ($m$)*], [*Intercept ($c$)*]),
    table.hline(),
    [Red], [0.9785], [-0.0095],
    [Near-Infrared (NIR)], [0.9548], [+0.0068],
    table.hline(),
  ),
  caption: [OLS Transformation Coefficients (Landsat 8 to 7 equivalent)],
)

*Harmonization Protocol:*
- *Landsat 5 & 7*: Treated as spectrally compatible (no transformation).
- *Landsat 8 & 9*: Transformed using $"L7"_("eq") = m dot "L8" + c$.

This correction reduces cross-sensor discontinuities that could otherwise be misinterpreted as vegetation trends.

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
The USGS explicitly states that "dense vegetation such as that found in temperate and tropical forests... typically exhibits NDVI values of approximately 0.6 to 0.9" (USGS, 2025). This threshold is the primary basis for defining the "Dense Canopy" class.

=== Sparse Vegetation (0.2 – 0.4)
USGS defines "shrub and grassland" as typically falling between 0.2 and 0.5 (USGS, 2025). This methodology uses 0.2–0.4 for sparse and 0.4–0.6 for transitional to segment this broad range.

= Trend Analysis

== Linear Trend Computation

Trends are calculated using Ordinary Least Squares (OLS) regression on annual summer median composites.

$ beta = (n sum x y - sum x sum y) / (n sum x^2 - (sum x)^2) $

Where $beta$ is the slope (change in NDVI per year) and $n$ is 40 years.

== Trend Significance Thresholds

To separate natural variability from significant change, a slope threshold of ±0.005 NDVI/year is applied.

#figure(
  table(
    columns: (auto, auto),
    inset: 8pt,
    align: left,
    table.header([*Trend Class*], [*Slope Threshold*]),
    table.hline(),
    [Gaining], [> +0.005/yr],
    [Stable], [±0.005/yr],
    [Losing], [< -0.005/yr],
    table.hline(),
  ),
  caption: [Trend Significance Thresholds],
)

=== Scientific Basis
This threshold is derived from Peng & Gong (2025), whose analysis of spatiotemporal NDVI changes classified slopes between 0.005 and 0.016 as "moderate improvement" and defined the stable range as -0.007 to 0.005. This provides a peer-reviewed basis for the cutoff.

= Classification Taxonomy

The system intersects absolute state (NDVI) with directional trend (Slope) to produce 9 mutually exclusive classes:

== Trend-Driven Classes (> ±0.005/yr)

#figure(
  table(
    columns: (auto, auto, auto),
    inset: 8pt,
    align: left,
    table.header([*Class*], [*Definition*], [*Transition Logic*]),
    table.hline(),
    [Canopy Establishment], [Sparse/Bare → Dense], [Crossed 0.6 threshold],
    [Emerging Biomass], [Sparse → Trans (+Gain)], [Early recovery phase],
    [Canopy Thickening], [Trans → Dense (+Gain)], [Maturation phase],
    [Canopy Densification], [Dense → Dense (+Gain)], [Increasing biomass],
    [Canopy Loss], [Dense → Sparse/Bare], [State collapse],
    [Canopy Thinning], [Dense → Trans (-Loss)], [Gradual degradation],
    table.hline(),
  ),
  caption: [Trend-Driven Classification Matrix],
)

== Stable/Edge Classes (< ±0.005/yr)

Captures slow transitions often found at forest ecotones:
- *Edge Expansion*: Sparse → Trans (Stable)
- *Edge Colonization*: Trans → Dense (Stable)
- *Edge Retreat*: Dense → Trans (Stable)

= Limitations & Caveats

Users must acknowledge the following limitations when interpreting results:

1.  *Threshold Universality*: While the 0.6 threshold is standard for temperate/tropical forests (USGS, 2025), boreal or dryland forests may require lower thresholds (e.g., 0.5).
2.  *Linearity Assumption*: Vegetation recovery often follows a sigmoid curve. Linear regression provides an average rate but may underestimate rapid recovery phases.
3.  *Sensor Homogeneity*: Despite harmonization (Roy et al., 2016), minor spectral differences between Landsat generations may influence trend calculations in subtle ways.

= Code Availability

The complete source code, including the Google Earth Engine script and documentation, is available at:

#link("https://github.com/gbrlpzz/ndvi-vegetation-cover-change")

= References

#set par(hanging-indent: 2em)

Peng, Y., & Gong, H. (2025). Analysis of Spatiotemporal Changes in NDVI-Derived Vegetation Index and Its Influencing Factors in Kunming City (2000 to 2020). _Forests_, 16(12), 1781. https://doi.org/10.3390/f16121781

Pettorelli, N., Vik, J. O., Mysterud, A., Gaillard, J. M., Tucker, C. J., & Stenseth, N. C. (2005). Using the satellite-derived NDVI to assess ecological responses to environmental change. _Trends in Ecology & Evolution_, 20(9), 503–510. https://doi.org/10.1016/j.tree.2005.05.011

Roy, D. P., Kovalskyy, V., Zhang, H. K., Vermote, E. F., Yan, L., Kumar, S. S., & Egorov, A. (2016). Characterization of Landsat-7 to Landsat-8 reflective wavelength and normalized difference vegetation index continuity. _Remote Sensing of Environment_, 185, 57–70. https://doi.org/10.1016/j.rse.2015.12.024

U.S. Geological Survey (USGS). (2025). _NDVI, the Foundation for Remote Sensing Phenology_. Retrieved from https://www.usgs.gov/special-topics/remote-sensing-phenology/science/ndvi-foundation-remote-sensing-phenology

#v(1cm)
#align(center)[
  #text(size: 9pt, fill: gray)[
    Document version 2.3 | December 2025 | © 2025 Gabriele Pizzi
  ]
]
