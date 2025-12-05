// ============================================================================
// VEGETATION COVER CHANGE DETECTION VIA NDVI TREND ANALYSIS ON GEE v2.0.0
// Multi-Decadal Landsat Collection 2 Analysis (Inter-calibrated)
//
// Author: Gabriele Pizzi (2025)
// License: Apache-2.0
// DOI: 10.5281/zenodo.17831391
// Repository: https://github.com/gbrlpzz/ndvi-vegetation-cover-change
// ============================================================================
//
// CLASSIFICATION TAXONOMY:
// -------------------------
// NDVI Thresholds:     Dense ≥0.6* | Trans 0.4-0.6* | Sparse 0.2-0.4* | Bare <0.2*
//                      (*Adjustable via SENSITIVITY_ADJUSTMENT)
// Trend Thresholds:    Gaining >+0.005/yr | Losing <-0.005/yr | Stable ±0.005
//                      (All trends filtered by Mann-Kendall Significance p < 0.05)
//
// CHANGE CLASSES (State-Driven):
//   1. Canopy Loss           - Dense → Sparse/Bare
//   2. Degradation           - Dense → Trans
//   3. Emerging Biomass      - Sparse → Trans
//   4. Maturation            - Trans → Dense
//   5. Canopy Densification  - Dense → Dense (gaining trend)
//   6. Canopy Establishment  - Sparse → Dense (epoch tracked)
//   7. Sparse Accumulation   - Sparse → Sparse (gaining trend)
//   8. Trans. Accumulation   - Trans → Trans (gaining trend)
//
// See /docs/methodology.pdf for full documentation.
// ============================================================================

// 1. CONFIGURATION

// Time Series Parameters
var startYear = 1985;
var endYear = 2025; // Update this to extend analysis
var recentYearStart = endYear - 10; // Dynamic 10-year trend comparison window

// SENSITIVITY ANALYSIS
// Set to a value (e.g. ±0.05) to test threshold stability. Default 0.
var SENSITIVITY_ADJUSTMENT = 0.0;

// Seasonality Parameters (Months 1-12)
// Default: June (6) to Sept (9) (Seasonal Window - Configurable)
var START_MONTH = 6;
var END_MONTH = 9;

// NDVI Thresholds (Literature-validated + Sensitivity)
var DENSE_CANOPY = 0.6 + SENSITIVITY_ADJUSTMENT;      // Dense forest canopy
var TRANSITIONAL = 0.4 + SENSITIVITY_ADJUSTMENT;      // Transitional woodland-shrub
var SPARSE = 0.2 + SENSITIVITY_ADJUSTMENT;            // Sparse vegetation / open land

// Trend Thresholds (Source: Peng & Gong, 2025)
var GAINING_SLOPE = 0.005;   // Active biomass accumulation
var LOSING_SLOPE = -0.005;   // Active biomass decline

// Region of Interest
// Region of Interest
var roi = roi || null;
if (typeof roi === 'undefined' || roi === null) {
  // Default fallback if undefined
  // Fallback: A known deforestation hotspot in Rondonia, Brazil
  roi = ee.Geometry.Polygon([
    [-63.0, -10.0], [-62.5, -10.0], [-62.5, -9.5], [-63.0, -9.5]
  ]);
  print('⚠️ WARNING: No ROI defined. Using default reproducible area (Rondonia, Brazil).');
  Map.centerObject(roi, 10);
} else {
  Map.centerObject(roi);
}

// 2. DATA PROCESSING

function maskL57(image) {
  var qa = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 3).eq(0).and(qa.bitwiseAnd(1 << 4).eq(0));
  return image.updateMask(mask)
    .select(['SR_B3', 'SR_B4'], ['Red', 'NIR'])
    .multiply(0.0000275).add(-0.2)
    .set('system:time_start', image.get('system:time_start'));
}

function maskL89(image) {
  var qa = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 3).eq(0).and(qa.bitwiseAnd(1 << 4).eq(0));
  // USGS Collection 2 Level-2 Surface Reflectance is already inter-calibrated
  // We strictly select and rename bands to ensure compatibility with L5/7 bands
  return image.select(['SR_B4', 'SR_B5'], ['Red', 'NIR'])
    .multiply(0.0000275).add(-0.2) // Scale to reflectance
    .updateMask(mask)
    .set('system:time_start', image.get('system:time_start'));
}

var l5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2").filterBounds(roi).map(maskL57);
var l7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2").filterBounds(roi).map(maskL57);
var l8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2").filterBounds(roi).map(maskL89);
var l9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2").filterBounds(roi).map(maskL89);
var fullCollection = l5.merge(l7).merge(l8).merge(l9);

// 3. COMPUTE NDVI STATES

// 94. COMPUTE NDVI STATES
// 95. 
function getSeasonalComposite(startDate, endDate) {
  return fullCollection.filterDate(startDate, endDate)
    .filter(ee.Filter.calendarRange(START_MONTH, END_MONTH, 'month'))
    .median()
    .normalizedDifference(['NIR', 'Red']);
}

// Baseline (startYear to startYear+4) and Current (endYear-4 to endYear) states
// Using 5-year averages for stable state definition
// Baseline (startYear to startYear+4) and Current (endYear-4 to endYear) states
// Using 5-year averages for stable state definition
var startNDVI = getSeasonalComposite(startYear + '-01-01', (startYear + 4) + '-12-31');
var endNDVI = getSeasonalComposite((endYear - 4) + '-01-01', endYear + '-12-31');

// Vegetation class: 1=Dense, 2=Transitional, 3=Sparse, 4=Bare
function classifyNDVI(ndvi) {
  return ee.Image(4)
    .where(ndvi.gte(SPARSE), 3)
    .where(ndvi.gte(TRANSITIONAL), 2)
    .where(ndvi.gte(DENSE_CANOPY), 1);
}

var startClass = classifyNDVI(startNDVI);
var endClass = classifyNDVI(endNDVI);

// 4. TREND ANALYSIS (Full period)

// 4. TREND ANALYSIS (Full period) with Statistical Significance

var trendCollection = fullCollection.filterDate(startYear + '-01-01', endYear + '-12-31')
  .filter(ee.Filter.calendarRange(START_MONTH, END_MONTH, 'month'))
  .map(function (img) {
    return img.normalizedDifference(['NIR', 'Red']).rename('NDVI')
      .set('system:time_start', img.get('system:time_start'));
  });

// 1. Linear Slope (Magnitude of change)
var linearFit = trendCollection.select(['NDVI'])
  .map(function (img) {
    // Note: Divide by 31536000000 (365 days) is an approximation. 
    // For scientific precision with leap years, consider: img.date().get('year').add(img.date().getFraction('year'))
    var t = ee.Image.constant(img.get('system:time_start')).divide(31536000000).float();
    return img.addBands(t.rename('t'));
  })
  .select(['t', 'NDVI'])
  .reduce(ee.Reducer.linearFit());

var slope = linearFit.select('scale');

// 2. Statistical Significance (Kendall's Tau)
// Non-parametric test robust to noise and non-normality
var kendall = trendCollection.select('NDVI') // Reducer assumes time series input
  .reduce(ee.Reducer.kendallsCorrelation());

var pValue = kendall.select('NDVI_p-value');
var SIGNIFICANCE_LEVEL = 0.05; // 95% Confidence

// Recent Trend (Dynamic)
var recentTrendCollection = fullCollection.filterDate(recentYearStart + '-01-01', endYear + '-12-31')
  .filter(ee.Filter.calendarRange(START_MONTH, END_MONTH, 'month'))
  .map(function (img) {
    var ndvi = img.normalizedDifference(['NIR', 'Red']).rename('NDVI');
    var t = ee.Image.constant(img.get('system:time_start')).divide(31536000000).float().rename('t');
    return ndvi.addBands(t);
  });

var recentFit = recentTrendCollection.select(['t', 'NDVI']).reduce(ee.Reducer.linearFit());
var recentSlope = recentFit.select('scale').rename('recent_slope');

// Trend class: 1=Gaining, 2=Stable, 3=Losing
// MASKED by statistical significance (p < 0.05)
var significantTrend = pValue.lt(SIGNIFICANCE_LEVEL);

var trendClass = ee.Image(2) // Default to Stable
  .where(slope.gt(GAINING_SLOPE).and(significantTrend), 1)
  .where(slope.lt(LOSING_SLOPE).and(significantTrend), 3);

// 5. CHANGE CLASSIFICATION (Simplified Taxonomy)
// Primarily driven by State Change (Robust Median Comparison)
// Trend Analysis used only for intra-class dynamics (Densification)

// 1 = Canopy Loss        (Dense → Sparse/Bare)
// 2 = Degradation        (Dense → Transitional)
// 3 = Emerging Biomass   (Sparse → Transitional)
// 4 = Maturation         (Transitional → Dense)
// 5 = Canopy Densification (Dense → Dense + Gaining)
// 6 = Canopy Establishment (Sparse/Bare → Dense)
// 7 = Sparse Accumulation       (Sparse → Sparse + Gaining)
// 8 = Transitional Accumulation (Trans → Trans + Gaining)

var changeClass = ee.Image(0);

// 1. Canopy Loss
changeClass = changeClass.where(
  startClass.eq(1).and(endClass.gte(3)), 1);

// 2. Degradation (formerly Thinning/Edge Retreat)
changeClass = changeClass.where(
  startClass.eq(1).and(endClass.eq(2)), 2);

// 3. Emerging Biomass (formerly Emerging/Edge Expansion)
changeClass = changeClass.where(
  startClass.eq(3).and(endClass.eq(2)), 3);

// 4. Maturation (formerly Thickening/Edge Colonization)
changeClass = changeClass.where(
  startClass.eq(2).and(endClass.eq(1)), 4);

// 5. Canopy Densification
// Requires stable Dense state AND significant gaining trend
changeClass = changeClass.where(
  startClass.eq(1).and(endClass.eq(1)).and(trendClass.eq(1)), 5);

// 6. Canopy Establishment
// Strict definition: Sparse/Bare → Dense (New forest)
var establishmentMask = startClass.gte(3).and(endClass.eq(1));
changeClass = changeClass.where(establishmentMask, 6);

// 7. Sparse Accumulation (Sparse → Sparse + Gaining)
var sparseStable = startClass.eq(3).and(endClass.eq(3));
var isGaining = trendClass.eq(1).or(recentSlope.gt(GAINING_SLOPE));
changeClass = changeClass.where(sparseStable.and(isGaining), 7);

// 8. Transitional Accumulation (Trans → Trans + Gaining)
var transStable = startClass.eq(2).and(endClass.eq(2));
changeClass = changeClass.where(transStable.and(isGaining), 8);

changeClass = changeClass.rename('change_class');

// 6. CANOPY ESTABLISHMENT EPOCHS (Dynamic Generation)

var epochs = [];
// Generate 5-year epochs starting from startYear + 5 (since first 5 are baseline)
// Generate 5-year epochs starting from startYear + 5
// Logic: If the remaining tail is small (< 3 years), merge it into the final epoch
// to avoid orphan 1-year epochs (e.g., 2020-2025 implies a 6-year final bin).
for (var y = startYear + 5; y <= endYear;) {
  var nextStart = y + 5;
  var yearsRemaining = endYear - nextStart + 1;

  if (yearsRemaining < 3) {
    // Merge remainder into this epoch and finish
    epochs.push({ start: y, end: endYear, label: y });
    break;
  } else {
    // Standard 5-year bin
    epochs.push({ start: y, end: y + 4, label: y });
    y += 5;
  }
}

var epochCollection = ee.ImageCollection.fromImages(
  epochs.map(function (epoch) {
    var img = getSeasonalComposite(epoch.start + '-01-01', epoch.end + '-12-31');
    return img.gte(DENSE_CANOPY)
      .multiply(epoch.label)
      .selfMask()
      .toInt()
      .rename('epoch')
      .set('epoch', epoch.label);
  })
);

// Mask includes both Establishment (Sparse->Dense) and Maturation (Trans->Dense)
// We track the year for ANY pixel that ends up Dense but started non-Dense
var epochValidMask = startClass.gte(2).and(endClass.eq(1));
var establishmentEpoch = epochCollection.min().updateMask(epochValidMask);

// 7. TRAJECTORY PROJECTION (Sigmoid-based)
// For gaining areas, project years to reach dense canopy threshold
// Using linear extrapolation: years = (threshold - currentNDVI) / slope
// LIMITATION: Capped at 50 years to avoid unrealistic ecological recovery times


var yearsToThreshold = endNDVI.subtract(DENSE_CANOPY).abs()
  .divide(recentSlope.abs())
  .where(recentSlope.lte(0), 9999)  // No projection for non-gaining
  .where(endNDVI.gte(DENSE_CANOPY), 0)  // Already at threshold
  .clamp(0, 50)
  .rename('years_to_canopy');

// Only show for areas currently gaining and below threshold
var projectionMask = trendClass.eq(1).and(endClass.gt(1));
var yearsToCanopy = yearsToThreshold.updateMask(projectionMask);

// 8. VISUALIZATION

Map.centerObject(roi);

// Define names for 6 simplified classes
var classNames = {
  1: 'Canopy Loss',
  2: 'Degradation',
  3: 'Emerging Biomass',
  4: 'Maturation',
  5: 'Densification',
  6: 'Establishment',
  7: 'Sparse Accumulation',
  8: 'Transitional Accumulation'
};

var changeViz = {
  min: 1,
  max: 8,
  palette: [
    'D7191C', // 1. Loss (Deep Carmine)
    'FDAE61', // 2. Degradation (Soft Orange)
    'A6D96A', // 3. Emerging (Yellow-Green)
    '008837', // 4. Maturation (Solid Green)
    '00441B', // 5. Densification (Deep Green)
    '2C7BB6', // 6. Establishment (Medium Blue)
    'D9F0D3', // 7. Sparse Accumulation (Pale Green)
    '74C476'  // 8. Transitional Accumul. (Soft Green)
  ]
};
Map.addLayer(changeClass.updateMask(changeClass.gt(0)), changeViz, 'Vegetation Change');

// Establishment epochs (off by default)
// Dynamic visualization range based on generated epochs
var epochViz = {
  min: epochs.length > 0 ? epochs[0].label : 1990,
  max: epochs.length > 0 ? epochs[epochs.length - 1].label : 2020,
  palette: ['08306b', '2171b5', '4eb3d3', '7fcdbb', 'c7e9b4', 'ffffb2', 'fd8d3c']
};
Map.addLayer(establishmentEpoch, epochViz, 'Canopy Gain Epoch (Est. + Mat.)', false);

// Years to canopy projection (off by default)
var projViz = {
  min: 0,
  max: 30,
  palette: ['00FF00', 'FFFF00', 'FF0000'] // Green (soon) to Red (distant)
};
Map.addLayer(yearsToCanopy, projViz, 'Years to Dense Canopy (Theoretical)', false);

// 8b. SUPPLEMENTARY DATA LAYERS (User Requested)

// Statistical Trend Class (All vegetation types)
// Shows significant trends (p < 0.05) even in Sparse/Transitional areas
var trendViz = {
  min: 1,
  max: 3,
  palette: ['228B22', 'CCCCCC', 'FF0000'] // 1:Gain(Green), 2:Stable(Grey), 3:Loss(Red)
};
// Mask out non-significant pixels for cleaner visualization
var significantTrendLayer = trendClass.updateMask(significantTrend);
Map.addLayer(significantTrendLayer, trendViz, 'Statistical Trends (All Classes)', false);

// Trend Magnitude (Slope)
var slopeViz = {
  min: -0.015,
  max: 0.015,
  palette: ['a50026', 'd73027', 'f46d43', 'fdae61', 'fee08b', 'ffffbf', 'd9ef8b', 'a6d96a', '66bd63', '1a9850', '006837']
};
Map.addLayer(slope.updateMask(significantTrend), slopeViz, 'Trend Magnitude (Slope)', false);

// DYNAMIC LEGEND - Updates based on visible layer
var legend = ui.Panel({
  style: { position: 'bottom-left', padding: '10px 14px', backgroundColor: 'white', maxWidth: '280px' }
});
Map.add(legend);

function makeRow(color, name, desc) {
  var colorBox = ui.Label({ style: { backgroundColor: '#' + color, padding: '10px', margin: '0 8px 6px 0' } });
  var textPanel = ui.Panel({
    widgets: [
      ui.Label({ value: name, style: { fontWeight: 'bold', fontSize: '11px', margin: '0' } }),
      ui.Label({ value: desc, style: { fontSize: '9px', color: '555555', margin: '0' } })
    ],
    style: { margin: '0' }
  });
  return ui.Panel({ widgets: [colorBox, textPanel], layout: ui.Panel.Layout.Flow('horizontal'), style: { margin: '0 0 2px 0' } });
}

function updateLegend(layerName) {
  legend.clear();

  if (layerName === 'Canopy Gain Epoch (Est. + Mat.)') {
    // Epoch legend
    legend.add(ui.Label({ value: 'Canopy Gain Epoch', style: { fontWeight: 'bold', fontSize: '14px', margin: '0 0 4px 0' } }));
    legend.add(ui.Label({ value: 'When area first reached dense canopy', style: { fontSize: '10px', color: '666666', margin: '0 0 10px 0' } }));

    // Dynamically generate legend rows from epochs array
    var palette = epochViz.palette;
    for (var i = 0; i < epochs.length; i++) {
      var e = epochs[i];
      // Cycle through palette if more epochs than colors
      var color = palette[i % palette.length];
      var label = e.start + '-' + e.end;
      var desc = (i === 0) ? 'Earliest' : ((i === epochs.length - 1) ? 'Latest' : '');
      legend.add(makeRow(color, label, desc));
    }

  } else if (layerName === 'Years to Dense Canopy (Theoretical)') {
    // Projection legend
    legend.add(ui.Label({ value: 'Years to Dense Canopy', style: { fontWeight: 'bold', fontSize: '14px', margin: '0 0 4px 0' } }));
    legend.add(ui.Label({ value: 'Theoretical Linear Projection', style: { fontSize: '10px', color: '666666', margin: '0 0 10px 0' } }));
    legend.add(makeRow('00FF00', '0-5 years', 'Imminent'));
    legend.add(makeRow('7FFF00', '5-10 years', ''));
    legend.add(makeRow('FFFF00', '10-20 years', ''));
    legend.add(makeRow('FF7F00', '20-30 years', ''));
    legend.add(makeRow('FF7F00', '20-30 years', ''));
    legend.add(makeRow('FF0000', '>30 years', 'Distant'));

  } else if (layerName === 'Statistical Trends (All Classes)') {
    // Trend legend
    legend.add(ui.Label({ value: 'Statistical Trends (MK Test)', style: { fontWeight: 'bold', fontSize: '14px', margin: '0 0 4px 0' } }));
    legend.add(ui.Label({ value: 'Significant trends (p < 0.05) across all strata', style: { fontSize: '10px', color: '666666', margin: '0 0 10px 0' } }));
    legend.add(makeRow('228B22', 'Gaining', '> +0.005 NDVI/yr'));
    legend.add(makeRow('FF0000', 'Losing', '< -0.005 NDVI/yr'));
    legend.add(makeRow('CCCCCC', 'Stable', 'No sig. trend or low slope'));

  } else if (layerName === 'Trend Magnitude (Slope)') {
    legend.add(ui.Label({ value: 'Trend Magnitude', style: { fontWeight: 'bold', fontSize: '14px', margin: '0 0 4px 0' } }));
    legend.add(ui.Label({ value: 'Slope of significant linear trends', style: { fontSize: '10px', color: '666666', margin: '0 0 10px 0' } }));
    // Simple gradient representation
    var gradient = ui.Panel({
      style: {
        width: '200px', height: '20px', margin: '0 0 10px 0',
        backgroundImage: 'linear-gradient(to right, #a50026, #ffffbf, #006837)'
      }
    });
    var labels = ui.Panel({
      widgets: [ui.Label('-0.015', { margin: '0', fontSize: '10px' }), ui.Label('0', { margin: '0 auto', fontSize: '10px' }), ui.Label('+0.015', { margin: '0', fontSize: '10px' })],
      layout: ui.Panel.Layout.Flow('horizontal'), style: { width: '200px' }
    });
    legend.add(gradient);
    legend.add(labels);

  } else {
    // Default: Vegetation Change legend
    legend.add(ui.Label({ value: 'Vegetation Cover Change', style: { fontWeight: 'bold', fontSize: '14px', margin: '0 0 4px 0' } }));
    var duration = endYear - startYear;
    legend.add(ui.Label({ value: duration + '-Year Analysis (' + startYear + '-' + endYear + ')', style: { fontSize: '10px', color: '666666', margin: '0 0 10px 0' } }));

    legend.add(ui.Label({ value: 'State Transitions (Robust)', style: { fontSize: '9px', color: '666666', margin: '0 0 4px 0' } }));
    legend.add(makeRow('D7191C', 'Canopy Loss', 'Dense → Sparse/Bare'));
    legend.add(makeRow('FDAE61', 'Degradation', 'Dense → Transitional'));
    legend.add(makeRow('A6D96A', 'Emerging Biomass', 'Sparse → Transitional'));
    legend.add(makeRow('008837', 'Maturation', 'Transitional → Dense'));
    legend.add(makeRow('2C7BB6', 'Establishment', 'Sparse → Dense'));

    legend.add(ui.Label({ value: 'Trend Confirmed', style: { fontSize: '9px', color: '666666', margin: '8px 0 4px 0' } }));
    legend.add(makeRow('00441B', 'Densification', 'Dense → Dense (+Gain)'));
    legend.add(makeRow('74C476', 'Trans. Accumulation', 'Trans → Trans (+Gain)'));
    legend.add(makeRow('D9F0D3', 'Sparse Accumulation', 'Sparse → Sparse (+Gain)'));

    var footerPanel = ui.Panel({ style: { margin: '10px 0 0 0', padding: '8px 0 0 0' } });
    footerPanel.add(ui.Label({ value: 'Dense ≥' + DENSE_CANOPY + ' | Trans ' + TRANSITIONAL + '-' + DENSE_CANOPY + ' | Sparse ' + SPARSE + '-' + TRANSITIONAL, style: { fontSize: '9px', color: '888888' } }));
    legend.add(footerPanel);
  }
}

// Initialize with default legend
updateLegend('Vegetation Change');

// Note: Legend updates automatically when you click map (checks visible layer)

// Create inspector panel (floating overlay like legend)
var inspectorPanel = ui.Panel({
  style: {
    position: 'bottom-right',
    padding: '8px',
    backgroundColor: 'white',
    width: '350px',
    shown: false
  }
});
Map.add(inspectorPanel);

// Lookup tables for labels
// classNames is already defined above
var vegNames = { 1: 'Dense Canopy', 2: 'Transitional', 3: 'Sparse', 4: 'Bare' };
var trendNames = { 1: 'Gaining', 2: 'Stable', 3: 'Losing' };

Map.style().set('cursor', 'crosshair');

// Inspector function
function updateInspector(coords) {
  var point = ee.Geometry.Point(coords.lon, coords.lat);

  // Update legend based on visible layer
  var layers = Map.layers();
  for (var i = 0; i < layers.length(); i++) {
    var layer = layers.get(i);
    if (layer.getShown()) {
      updateLegend(layer.getName());
      break;
    }
  }

  // Clear and show panel
  inspectorPanel.clear();
  inspectorPanel.style().set('shown', true);

  // Loading indicator
  inspectorPanel.add(ui.Label('Loading...', { fontStyle: 'italic' }));

  // Get values at point
  var values = ee.Image.cat([
    changeClass,
    startClass.rename('start_class'),
    endClass.rename('end_class'),
    trendClass.rename('trend_class'),
    slope.rename('slope'),
    pValue.rename('p_value'),
    recentSlope,
    endNDVI.rename('current_ndvi'),
    establishmentEpoch.rename('epoch'),
    yearsToCanopy.rename('years_proj')
  ]).reduceRegion({
    reducer: ee.Reducer.first(),
    geometry: point,
    scale: 30
  });

  // Create NDVI time series
  var ndviSeries = fullCollection
    .filter(ee.Filter.calendarRange(START_MONTH, END_MONTH, 'month'))
    .map(function (img) {
      var ndvi = img.normalizedDifference(['NIR', 'Red']);
      return ee.Feature(point, {
        'NDVI': ndvi.reduceRegion({
          reducer: ee.Reducer.first(),
          geometry: point,
          scale: 30
        }).get('nd'),
        'system:time_start': img.get('system:time_start')
      });
    });

  values.evaluate(function (res) {
    inspectorPanel.clear();

    if (!res) {
      inspectorPanel.add(ui.Label('No data available (masked or outside bounds)', { fontSize: '11px', color: 'gray' }));
      return;
    }

    // Header with close button
    var header = ui.Panel({
      widgets: [
        ui.Label('Point Analysis', { fontWeight: 'bold', fontSize: '14px' }),
        ui.Button({
          label: '✕', style: { padding: '0 6px' },
          onClick: function () { inspectorPanel.style().set('shown', false); }
        })
      ],
      layout: ui.Panel.Layout.Flow('horizontal'),
      style: { margin: '0 0 6px 0' }
    });
    inspectorPanel.add(header);

    // Coordinates
    inspectorPanel.add(ui.Label(
      'Location: ' + coords.lat.toFixed(5) + '°N, ' + coords.lon.toFixed(5) + '°E',
      { fontSize: '10px', color: '666666', margin: '0 0 8px 0' }
    ));

    // Current state section
    var currentYearStart = endYear - 4;
    inspectorPanel.add(ui.Label('Current State (' + currentYearStart + '-' + endYear + ')', { fontWeight: 'bold', fontSize: '11px', margin: '0 0 4px 0' }));

    var currentNDVI = res.current_ndvi || 0;
    var currentClass = vegNames[res.end_class] || 'Unknown';
    // Dynamic color coding based on active thresholds
    var ndviColor = currentNDVI >= DENSE_CANOPY ? '006400' : (currentNDVI >= TRANSITIONAL ? '90EE90' : (currentNDVI >= SPARSE ? 'AAAA00' : 'AA6600'));

    var currentPanel = ui.Panel({
      widgets: [
        ui.Label('NDVI: ' + currentNDVI.toFixed(3), { fontSize: '11px', fontWeight: 'bold' }),
        ui.Label(currentClass, { fontSize: '11px', backgroundColor: '#' + ndviColor, color: currentNDVI >= TRANSITIONAL ? 'white' : 'black', padding: '2px 6px' })
      ],
      layout: ui.Panel.Layout.Flow('horizontal'),
      style: { margin: '0 0 8px 0' }
    });
    inspectorPanel.add(currentPanel);

    // Trend section
    var slopeVal = res.slope || 0;
    var recentSlopeVal = res.recent_slope || 0;
    var slopeClass = slopeVal > GAINING_SLOPE ? 'Gaining' : (slopeVal < LOSING_SLOPE ? 'Losing' : 'Stable');
    var slopeColor = slopeVal > GAINING_SLOPE ? '228B22' : (slopeVal < LOSING_SLOPE ? 'CC0000' : '888888');
    var recentClass = recentSlopeVal > GAINING_SLOPE ? 'Gaining' : (recentSlopeVal < LOSING_SLOPE ? 'Losing' : 'Stable');
    var recentColor = recentSlopeVal > GAINING_SLOPE ? '228B22' : (recentSlopeVal < LOSING_SLOPE ? 'CC0000' : '888888');

    // Determine trend acceleration
    var trendAccel = '';
    var trendAccelColor = '888888';
    var ACCEL_THRESHOLD = 0.002;
    if (Math.abs(recentSlopeVal - slopeVal) < ACCEL_THRESHOLD) {
      trendAccel = '→ Consistent';
      trendAccelColor = '666666';
    } else if (recentSlopeVal > slopeVal + ACCEL_THRESHOLD) {
      trendAccel = '↑ Accelerating';
      trendAccelColor = '228B22';
    } else {
      trendAccel = '↓ Decelerating';
      trendAccelColor = 'CC6600';
    }

    inspectorPanel.add(ui.Label('Trend Analysis', { fontWeight: 'bold', fontSize: '11px', margin: '0 0 4px 0' }));

    // Long-term row
    var duration = endYear - startYear;
    var trend40Panel = ui.Panel({
      widgets: [
        ui.Label(duration + 'yr:', { fontSize: '10px', color: '666666' }),
        ui.Label((slopeVal * 1000).toFixed(2) + '×10⁻³/yr', { fontSize: '10px' }),
        ui.Label(slopeClass, { fontSize: '10px', color: slopeColor, fontWeight: 'bold' })
      ],
      layout: ui.Panel.Layout.Flow('horizontal'),
      style: { margin: '0 0 2px 0' }
    });
    inspectorPanel.add(trend40Panel);

    var pVal = res.p_value || 1.0;
    var sigLabel = pVal < 0.05 ? 'Significant (p<0.05)' : 'Not Significant (p≥0.05)';
    var sigColor = pVal < 0.05 ? '228B22' : '888888';

    inspectorPanel.add(ui.Label(sigLabel, { fontSize: '9px', color: sigColor, margin: '0 0 4px 0', fontStyle: 'italic' }));

    // Recent row
    var trend10Panel = ui.Panel({
      widgets: [
        ui.Label('10yr:', { fontSize: '10px', color: '666666' }),
        ui.Label((recentSlopeVal * 1000).toFixed(2) + '×10⁻³/yr', { fontSize: '10px' }),
        ui.Label(recentClass, { fontSize: '10px', color: recentColor, fontWeight: 'bold' })
      ],
      layout: ui.Panel.Layout.Flow('horizontal'),
      style: { margin: '0 0 2px 0' }
    });
    inspectorPanel.add(trend10Panel);

    // Trend Acceleration
    inspectorPanel.add(ui.Label(trendAccel, {
      fontSize: '10px', color: trendAccelColor, fontStyle: 'italic', margin: '0 0 8px 0'
    }));

    // State transition
    inspectorPanel.add(ui.Label('State Transition', { fontWeight: 'bold', fontSize: '11px', margin: '0 0 4px 0' }));
    inspectorPanel.add(ui.Label(
      startYear + ': ' + (vegNames[res.start_class] || '?') + '  →  ' + endYear + ': ' + (vegNames[res.end_class] || '?'),
      { fontSize: '11px', margin: '0 0 8px 0' }
    ));

    // Change classification
    inspectorPanel.add(ui.Label('Classification Result', { fontWeight: 'bold', fontSize: '11px', margin: '0 0 4px 0' }));
    if (res.change_class && res.change_class > 0) {
      var changeColors = ['D7191C', 'FDAE61', 'A6D96A', '008837', '00441B', '2C7BB6', 'D9F0D3', '74C476'];
      var changeColor = changeColors[res.change_class - 1] || 'CCCCCC';
      var lightTextClasses = [3, 7];
      var changePanel = ui.Panel({
        style: { backgroundColor: '#' + changeColor, padding: '6px 10px', margin: '0 0 4px 0' }
      });
      changePanel.add(ui.Label(classNames[res.change_class] || 'Unknown', {
        fontWeight: 'bold', fontSize: '12px',
        color: lightTextClasses.indexOf(res.change_class) >= 0 ? '000000' : 'FFFFFF',
        backgroundColor: '#' + changeColor
      }));
      inspectorPanel.add(changePanel);
    } else {
      inspectorPanel.add(ui.Label('No class change detected', { fontSize: '11px', fontStyle: 'italic', color: '888888' }));
    }

    // Canopy Status (unified epoch/projection)
    inspectorPanel.add(ui.Label('Dense Canopy Status', { fontWeight: 'bold', fontSize: '11px', margin: '8px 0 4px 0' }));

    var statusPanel = ui.Panel({ style: { margin: '0 0 4px 0' } });
    var startedDense = res.start_class === 1;
    var nowDense = res.current_ndvi >= DENSE_CANOPY;

    if (res.epoch) {
      // Dynamic epoch end calculation
      var epochEnd = res.epoch === epochs[epochs.length - 1].label ? endYear : res.epoch + 4;
      statusPanel.add(ui.Label('✓ Canopy Reached: ' + res.epoch + '-' + epochEnd, {
        fontSize: '11px', color: '228B22', fontWeight: 'bold'
      }));
      statusPanel.add(ui.Label('First reached NDVI ≥' + DENSE_CANOPY + ' (from non-dense)', {
        fontSize: '9px', color: '666666', fontStyle: 'italic'
      }));
    } else if (startedDense && nowDense) {
      statusPanel.add(ui.Label('✓ Established: pre-' + startYear, {
        fontSize: '11px', color: '228B22', fontWeight: 'bold'
      }));
      statusPanel.add(ui.Label('Dense canopy throughout study period', {
        fontSize: '9px', color: '666666', fontStyle: 'italic'
      }));
    } else if (!startedDense && nowDense && slopeVal > 0) {
      // Transitioned TO Dense
      var yearsAgo = (currentNDVI - DENSE_CANOPY) / slopeVal;
      var crossYear = Math.round(endYear - yearsAgo);
      statusPanel.add(ui.Label('✓ Established: ~' + crossYear, {
        fontSize: '11px', color: '228B22', fontWeight: 'bold'
      }));
      statusPanel.add(ui.Label('Estimated crossing from ' + (vegNames[res.start_class] || 'transitional'), {
        fontSize: '9px', color: '666666', fontStyle: 'italic'
      }));
    } else if (!startedDense && nowDense) {
      statusPanel.add(ui.Label('✓ Established: ~2020s', {
        fontSize: '11px', color: '228B22', fontWeight: 'bold'
      }));
      statusPanel.add(ui.Label('Recently reached threshold', {
        fontSize: '9px', color: '666666', fontStyle: 'italic'
      }));
    } else if (slopeVal > 0) {
      // Gaining
      var yearsNeeded = (DENSE_CANOPY - currentNDVI) / recentSlopeVal;
      if (yearsNeeded <= 0) {
        statusPanel.add(ui.Label('✓ At threshold', {
          fontSize: '11px', color: '228B22', fontWeight: 'bold'
        }));
      } else if (yearsNeeded <= 50) {
        var projYear = endYear + Math.round(yearsNeeded);
        statusPanel.add(ui.Label('↗ Projected: ~' + projYear + ' (Theoretical)', {
          fontSize: '11px', color: '0066CC', fontWeight: 'bold'
        }));
        statusPanel.add(ui.Label('Est. ' + Math.round(yearsNeeded) + ' years (Linear Model)', {
          fontSize: '9px', color: '666666', fontStyle: 'italic'
        }));
      } else {
        statusPanel.add(ui.Label('↗ Projected: ~' + (endYear + Math.round(yearsNeeded)) + ' (slow)', {
          fontSize: '11px', color: '888888'
        }));
        statusPanel.add(ui.Label('~' + Math.round(yearsNeeded) + ' years at current rate', {
          fontSize: '9px', color: '999999', fontStyle: 'italic'
        }));
      }
    } else if (slopeVal < LOSING_SLOPE) {
      statusPanel.add(ui.Label('↘ Declining (no projection)', {
        fontSize: '11px', color: 'CC0000'
      }));
      statusPanel.add(ui.Label('Strong negative trend', {
        fontSize: '9px', color: '999999', fontStyle: 'italic'
      }));
    } else if (slopeVal < 0) {
      statusPanel.add(ui.Label('↘ Slight decline', {
        fontSize: '11px', color: 'CC6600'
      }));
      statusPanel.add(ui.Label('Trend negative but below ' + LOSING_SLOPE + '/yr', {
        fontSize: '9px', color: '999999', fontStyle: 'italic'
      }));
    } else {
      statusPanel.add(ui.Label('— Stable (no change)', {
        fontSize: '11px', color: '888888'
      }));
      statusPanel.add(ui.Label('No measurable NDVI trend', {
        fontSize: '9px', color: '999999', fontStyle: 'italic'
      }));
    }
    inspectorPanel.add(statusPanel);

    // NDVI Chart
    var chart = ui.Chart.feature.byFeature({
      features: ndviSeries,
      xProperty: 'system:time_start',
      yProperties: ['NDVI']
    }).setChartType('ScatterChart')
      .setOptions({
        title: 'NDVI Trend (' + startYear + '-' + endYear + ')',
        titleTextStyle: { fontSize: 11, bold: true },
        hAxis: { title: '', format: 'yyyy', textStyle: { fontSize: 9 } },
        vAxis: {
          title: '',
          viewWindow: { min: 0, max: 1 },
          textStyle: { fontSize: 9 },
          gridlines: { count: 3 }
        },
        pointSize: 2,
        legend: { position: 'none' },
        series: { 0: { color: '333333' } },
        trendlines: { 0: { type: 'linear', color: 'FF0000', lineWidth: 2, showR2: true } },
        chartArea: { width: '85%', height: '70%' },
        height: 180
      });

    inspectorPanel.add(chart);

    // Threshold reference
    inspectorPanel.add(ui.Label(
      'Thresholds: Dense≥' + DENSE_CANOPY + ' | Trans ' + TRANSITIONAL + '-' + DENSE_CANOPY + ' | Sparse ' + SPARSE + '-' + TRANSITIONAL,
      { fontSize: '9px', color: '888888', margin: '4px 0 0 0' }
    ));

    // Sensitivity Warning
    if (SENSITIVITY_ADJUSTMENT !== 0) {
      inspectorPanel.add(ui.Label(
        '⚠️ Sensitivity Analysis Active (Adj: ' + (SENSITIVITY_ADJUSTMENT > 0 ? '+' : '') + SENSITIVITY_ADJUSTMENT + ')',
        { fontSize: '9px', color: 'BC8F8F', margin: '2px 0 0 0', fontWeight: 'bold' }
      ));
    }
  });
}

Map.onClick(updateInspector);

// 10. EXPORT

Export.image.toDrive({
  image: changeClass.byte(),
  description: 'Export_Change_Classes_' + startYear + '_' + endYear,
  scale: 30,
  region: roi,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

Export.image.toDrive({
  image: establishmentEpoch.unmask(0).short(),
  description: 'Export_Establishment_Epoch_' + startYear + '_' + endYear,
  scale: 30,
  region: roi,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

Export.image.toDrive({
  image: slope.float(),
  description: 'Export_Trend_Slope_' + startYear + '_' + endYear,
  scale: 30,
  region: roi,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

Export.image.toDrive({
  image: trendClass.byte(),
  description: 'Export_Statistical_Trend_Class_' + startYear + '_' + endYear,
  scale: 30,
  region: roi,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});
