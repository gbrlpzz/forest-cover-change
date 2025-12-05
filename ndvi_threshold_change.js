// ============================================================================
// VEGETATION COVER CHANGE DETECTION v2.0
// 40-Year Landsat Analysis (1985-2025)
// ============================================================================
//
// CLASSIFICATION TAXONOMY:
// -------------------------
// NDVI Thresholds:     Dense ≥0.6 | Trans 0.4-0.6 | Sparse 0.2-0.4 | Bare <0.2
// Trend Thresholds:    Gaining >+0.005/yr | Losing <-0.005/yr | Stable ±0.005
//
// CHANGE CLASSES (Strong Trend):
//   1. Canopy Loss        - Dense → Sparse/Bare
//   2. Canopy Thinning    - Dense → Trans (losing)
//   3. Emerging Biomass   - Sparse → Trans (gaining)
//   4. Canopy Thickening  - Trans → Dense (gaining)
//   5. Canopy Densification - Dense → Dense (gaining)
//   6. Canopy Establishment - Sparse → Dense (epoch tracked)
//
// EDGE CLASSES (Stable Trend):
//   7. Edge Expansion     - Sparse → Trans
//   8. Edge Colonization  - Trans → Dense
//   9. Edge Retreat       - Dense → Trans
//
// See /docs/methodology.md for full documentation.
// ============================================================================

// 1. CONFIGURATION

var startYear = 1985;
var endYear = 2025;

// NDVI Thresholds (Literature-validated)
// Source: Multiple peer-reviewed studies on vegetation classification
var DENSE_CANOPY = 0.6;      // Dense forest canopy (≥30% cover)
var TRANSITIONAL = 0.4;      // Transitional woodland-shrub
var SPARSE = 0.2;            // Sparse vegetation / open land

// Trend Thresholds (Source: MDPI Kunming study)
// Slope of 0.005 NDVI/year = significant change
var GAINING_SLOPE = 0.005;   // Active biomass accumulation
var LOSING_SLOPE = -0.005;   // Active biomass decline

// Region of Interest
var roi = roi || Map.getBounds(true);

// 2. DATA PROCESSING

function maskL57(image) {
  var qa = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 3).eq(0).and(qa.bitwiseAnd(1 << 4).eq(0));
  return image.updateMask(mask)
    .select(['SR_B4', 'SR_B3'], ['NIR', 'Red'])
    .multiply(0.0000275).add(-0.2)
    .set('system:time_start', image.get('system:time_start'));
}

function maskL89(image) {
  var qa = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 3).eq(0).and(qa.bitwiseAnd(1 << 4).eq(0));
  return image.updateMask(mask)
    .select(['SR_B5', 'SR_B4'], ['NIR', 'Red'])
    .multiply(0.0000275).add(-0.2)
    .set('system:time_start', image.get('system:time_start'));
}

var l5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2").filterBounds(roi).map(maskL57);
var l7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2").filterBounds(roi).map(maskL57);
var l8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2").filterBounds(roi).map(maskL89);
var l9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2").filterBounds(roi).map(maskL89);
var fullCollection = l5.merge(l7).merge(l8).merge(l9);

// 3. COMPUTE NDVI STATES

function getSummerNDVI(startDate, endDate) {
  return fullCollection.filterDate(startDate, endDate)
    .filter(ee.Filter.calendarRange(6, 9, 'month'))
    .median()
    .normalizedDifference(['NIR', 'Red']);
}

// Baseline (1985-1989) and Current (2021-2025) states
var startNDVI = getSummerNDVI('1985-01-01', '1989-12-31');
var endNDVI = getSummerNDVI('2021-01-01', '2025-12-31');

// Vegetation class: 1=Dense, 2=Transitional, 3=Sparse, 4=Bare
function classifyNDVI(ndvi) {
  return ee.Image(4)
    .where(ndvi.gte(SPARSE), 3)
    .where(ndvi.gte(TRANSITIONAL), 2)
    .where(ndvi.gte(DENSE_CANOPY), 1);
}

var startClass = classifyNDVI(startNDVI);
var endClass = classifyNDVI(endNDVI);

// 4. TREND ANALYSIS (Full 40-year period for robust trend)

var trendCollection = fullCollection.filterDate('1985-01-01', '2025-12-31')
  .filter(ee.Filter.calendarRange(6, 9, 'month'))
  .map(function (img) {
    var ndvi = img.normalizedDifference(['NIR', 'Red']).rename('NDVI');
    var t = ee.Image.constant(img.get('system:time_start')).divide(31536000000).float().rename('t');
    return ndvi.addBands(t);
  });

var linearFit = trendCollection.select(['t', 'NDVI']).reduce(ee.Reducer.linearFit());
var slope = linearFit.select('scale');
var intercept = linearFit.select('offset');

// Recent 10-year trend (2015-2025) for momentum comparison
var recentTrendCollection = fullCollection.filterDate('2015-01-01', '2025-12-31')
  .filter(ee.Filter.calendarRange(6, 9, 'month'))
  .map(function (img) {
    var ndvi = img.normalizedDifference(['NIR', 'Red']).rename('NDVI');
    var t = ee.Image.constant(img.get('system:time_start')).divide(31536000000).float().rename('t');
    return ndvi.addBands(t);
  });

var recentFit = recentTrendCollection.select(['t', 'NDVI']).reduce(ee.Reducer.linearFit());
var recentSlope = recentFit.select('scale').rename('recent_slope');

// Trend class: 1=Gaining, 2=Stable, 3=Losing
var trendClass = ee.Image(2)
  .where(slope.gt(GAINING_SLOPE), 1)
  .where(slope.lt(LOSING_SLOPE), 3);

// 5. CHANGE CLASSIFICATION (Neutral Terminology)
// 1 = Canopy Loss (Dense → Sparse/Bare)
// 2 = Canopy Thinning (Dense → Transitional + Losing)
// 3 = Emerging Biomass (Sparse → Transitional + Gaining)
// 4 = Canopy Thickening (Transitional → Dense + Gaining)
// 5 = Canopy Densification (Dense → Dense + Gaining)
// 6 = Canopy Establishment (Sparse/Bare → Dense, tracked by epoch)

var changeClass = ee.Image(0);

// Canopy Loss: Dense (1) → Sparse (3) or Bare (4)
changeClass = changeClass.where(
  startClass.eq(1).and(endClass.gte(3)), 1);

// Canopy Thinning: Dense (1) → Transitional (2) + Losing
changeClass = changeClass.where(
  startClass.eq(1).and(endClass.eq(2)).and(trendClass.eq(3)), 2);

// Emerging Biomass: Sparse (3) → Transitional (2) + Gaining
changeClass = changeClass.where(
  startClass.eq(3).and(endClass.eq(2)).and(trendClass.eq(1)), 3);

// Canopy Thickening: Transitional (2) → Dense (1) + Gaining
changeClass = changeClass.where(
  startClass.eq(2).and(endClass.eq(1)).and(trendClass.eq(1)), 4);

// Canopy Densification: Dense (1) → Dense (1) + Gaining
changeClass = changeClass.where(
  startClass.eq(1).and(endClass.eq(1)).and(trendClass.eq(1)), 5);

// Canopy Establishment: Sparse/Bare → Dense
var establishmentMask = startClass.gte(3).and(endClass.eq(1));
changeClass = changeClass.where(establishmentMask, 6);

// EDGE EXPANSION CLASSES (state transitions with stable trend)
// These capture gradual edge expansion without strong directional trend

// Edge Expansion: Sparse (3) → Transitional (2) + Stable (not already flagged as Emerging)
changeClass = changeClass.where(
  startClass.eq(3).and(endClass.eq(2)).and(changeClass.eq(0)), 7);

// Edge Colonization: Transitional (2) → Dense (1) + Stable (not already flagged)
changeClass = changeClass.where(
  startClass.eq(2).and(endClass.eq(1)).and(changeClass.eq(0)), 8);

// Edge Retreat: Dense (1) → Transitional (2) + Stable (not thinning with losing trend)
changeClass = changeClass.where(
  startClass.eq(1).and(endClass.eq(2)).and(changeClass.eq(0)), 9);

changeClass = changeClass.rename('change_class');

// 6. CANOPY ESTABLISHMENT EPOCHS

var epochs = [
  { start: 1990, end: 1994, label: 1990 },
  { start: 1995, end: 1999, label: 1995 },
  { start: 2000, end: 2004, label: 2000 },
  { start: 2005, end: 2009, label: 2005 },
  { start: 2010, end: 2014, label: 2010 },
  { start: 2015, end: 2019, label: 2015 },
  { start: 2020, end: 2025, label: 2020 }
];

var epochCollection = ee.ImageCollection.fromImages(
  epochs.map(function (epoch) {
    var img = getSummerNDVI(epoch.start + '-01-01', epoch.end + '-12-31');
    return img.gte(DENSE_CANOPY)
      .multiply(epoch.label)
      .selfMask()
      .toInt()
      .rename('epoch')
      .set('epoch', epoch.label);
  })
);

var establishmentEpoch = epochCollection.min().updateMask(establishmentMask);

// 7. TRAJECTORY PROJECTION (Sigmoid-based)
// For gaining areas, project years to reach dense canopy threshold
// Using linear extrapolation: years = (threshold - currentNDVI) / slope

var yearsToThreshold = endNDVI.subtract(DENSE_CANOPY).abs()
  .divide(slope.abs())
  .where(slope.lte(0), 9999)  // No projection for non-gaining
  .where(endNDVI.gte(DENSE_CANOPY), 0)  // Already at threshold
  .clamp(0, 50)
  .rename('years_to_canopy');

// Only show for areas currently gaining and below threshold
var projectionMask = trendClass.eq(1).and(endClass.gt(1));
var yearsToCanopy = yearsToThreshold.updateMask(projectionMask);

// 8. VISUALIZATION

Map.centerObject(roi);

// Main change layer
var changeViz = {
  min: 1,
  max: 9,
  palette: [
    'FF00FF', // 1: Canopy Loss (Magenta)
    'FFA500', // 2: Canopy Thinning (Orange)
    'ADFF2F', // 3: Emerging Biomass (Lime)
    '90EE90', // 4: Canopy Thickening (Light Green)
    '006400', // 5: Canopy Densification (Dark Green)
    '0000FF', // 6: Canopy Establishment (Blue)
    'FFFF00', // 7: Edge Expansion (Yellow)
    '00CED1', // 8: Edge Colonization (Turquoise)
    'DDA0DD'  // 9: Edge Retreat (Plum)
  ]
};
Map.addLayer(changeClass.updateMask(changeClass.gt(0)), changeViz, 'Vegetation Change');

// Establishment epochs (off by default)
var epochViz = {
  min: 1990,
  max: 2020,
  palette: ['08306b', '2171b5', '4eb3d3', '7fcdbb', 'c7e9b4', 'ffffb2', 'fd8d3c']
};
Map.addLayer(establishmentEpoch, epochViz, 'Canopy Establishment Epoch', false);

// Years to canopy projection (off by default)
var projViz = {
  min: 0,
  max: 30,
  palette: ['00FF00', 'FFFF00', 'FF0000'] // Green (soon) to Red (distant)
};
Map.addLayer(yearsToCanopy, projViz, 'Years to Dense Canopy (Projection)', false);

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

function makeGradient(colors, labels) {
  var gradientPanel = ui.Panel({ layout: ui.Panel.Layout.Flow('horizontal'), style: { margin: '4px 0' } });
  for (var i = 0; i < colors.length; i++) {
    gradientPanel.add(ui.Label({
      value: labels[i],
      style: { backgroundColor: '#' + colors[i], color: i < colors.length / 2 ? 'white' : 'black', padding: '4px 8px', fontSize: '9px', margin: '0' }
    }));
  }
  return gradientPanel;
}

function updateLegend(layerName) {
  legend.clear();

  if (layerName === 'Canopy Establishment Epoch') {
    // Epoch legend
    legend.add(ui.Label({ value: 'Canopy Establishment Epoch', style: { fontWeight: 'bold', fontSize: '14px', margin: '0 0 4px 0' } }));
    legend.add(ui.Label({ value: 'When area first reached dense canopy', style: { fontSize: '10px', color: '666666', margin: '0 0 10px 0' } }));
    legend.add(makeRow('08306b', '1990-1994', 'Earliest'));
    legend.add(makeRow('2171b5', '1995-1999', ''));
    legend.add(makeRow('4eb3d3', '2000-2004', ''));
    legend.add(makeRow('7fcdbb', '2005-2009', ''));
    legend.add(makeRow('c7e9b4', '2010-2014', ''));
    legend.add(makeRow('ffffb2', '2015-2019', ''));
    legend.add(makeRow('fd8d3c', '2020-2025', 'Latest'));

  } else if (layerName === 'Years to Dense Canopy (Projection)') {
    // Projection legend
    legend.add(ui.Label({ value: 'Years to Dense Canopy', style: { fontWeight: 'bold', fontSize: '14px', margin: '0 0 4px 0' } }));
    legend.add(ui.Label({ value: 'Linear projection for gaining areas', style: { fontSize: '10px', color: '666666', margin: '0 0 10px 0' } }));
    legend.add(makeRow('00FF00', '0-5 years', 'Imminent'));
    legend.add(makeRow('7FFF00', '5-10 years', ''));
    legend.add(makeRow('FFFF00', '10-20 years', ''));
    legend.add(makeRow('FF7F00', '20-30 years', ''));
    legend.add(makeRow('FF0000', '>30 years', 'Distant'));

  } else {
    // Default: Vegetation Change legend
    legend.add(ui.Label({ value: 'Vegetation Cover Change', style: { fontWeight: 'bold', fontSize: '14px', margin: '0 0 4px 0' } }));
    legend.add(ui.Label({ value: '40-Year Analysis (1985-2025)', style: { fontSize: '10px', color: '666666', margin: '0 0 10px 0' } }));

    legend.add(ui.Label({ value: 'With Strong Trend (>±0.005/yr)', style: { fontSize: '9px', color: '666666', margin: '0 0 4px 0' } }));
    legend.add(makeRow('FF00FF', 'Canopy Loss', 'Dense → Sparse/Bare'));
    legend.add(makeRow('FFA500', 'Canopy Thinning', 'Dense → Trans (losing)'));
    legend.add(makeRow('ADFF2F', 'Emerging Biomass', 'Sparse → Trans (gaining)'));
    legend.add(makeRow('90EE90', 'Canopy Thickening', 'Trans → Dense (gaining)'));
    legend.add(makeRow('006400', 'Canopy Densification', 'Dense → Dense (gaining)'));
    legend.add(makeRow('0000FF', 'Canopy Establishment', 'Sparse → Dense'));

    legend.add(ui.Label({ value: 'With Stable Trend (<±0.005/yr)', style: { fontSize: '9px', color: '666666', margin: '8px 0 4px 0' } }));
    legend.add(makeRow('FFFF00', 'Edge Expansion', 'Sparse → Trans'));
    legend.add(makeRow('00CED1', 'Edge Colonization', 'Trans → Dense'));
    legend.add(makeRow('DDA0DD', 'Edge Retreat', 'Dense → Trans'));

    var footerPanel = ui.Panel({ style: { margin: '10px 0 0 0', padding: '8px 0 0 0' } });
    footerPanel.add(ui.Label({ value: 'Dense ≥0.6 | Trans 0.4-0.6 | Sparse 0.2-0.4', style: { fontSize: '9px', color: '888888' } }));
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
var classNames = {
  1: 'Canopy Loss',
  2: 'Canopy Thinning',
  3: 'Emerging Biomass',
  4: 'Canopy Thickening',
  5: 'Canopy Densification',
  6: 'Canopy Establishment',
  7: 'Edge Expansion',
  8: 'Edge Colonization',
  9: 'Edge Retreat'
};
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
    .filter(ee.Filter.calendarRange(6, 9, 'month'))
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
    inspectorPanel.add(ui.Label('Current State (2021-2025)', { fontWeight: 'bold', fontSize: '11px', margin: '0 0 4px 0' }));

    var currentNDVI = res.current_ndvi || 0;
    var currentClass = vegNames[res.end_class] || 'Unknown';
    var ndviColor = currentNDVI >= 0.6 ? '006400' : (currentNDVI >= 0.4 ? '90EE90' : (currentNDVI >= 0.2 ? 'AAAA00' : 'AA6600'));

    var currentPanel = ui.Panel({
      widgets: [
        ui.Label('NDVI: ' + currentNDVI.toFixed(3), { fontSize: '11px', fontWeight: 'bold' }),
        ui.Label(currentClass, { fontSize: '11px', backgroundColor: '#' + ndviColor, color: currentNDVI >= 0.4 ? 'white' : 'black', padding: '2px 6px' })
      ],
      layout: ui.Panel.Layout.Flow('horizontal'),
      style: { margin: '0 0 8px 0' }
    });
    inspectorPanel.add(currentPanel);

    // Trend section
    var slopeVal = res.slope || 0;
    var recentSlopeVal = res.recent_slope || 0;
    var slopeClass = slopeVal > 0.005 ? 'Gaining' : (slopeVal < -0.005 ? 'Losing' : 'Stable');
    var slopeColor = slopeVal > 0.005 ? '228B22' : (slopeVal < -0.005 ? 'CC0000' : '888888');
    var recentClass = recentSlopeVal > 0.005 ? 'Gaining' : (recentSlopeVal < -0.005 ? 'Losing' : 'Stable');
    var recentColor = recentSlopeVal > 0.005 ? '228B22' : (recentSlopeVal < -0.005 ? 'CC0000' : '888888');

    // Determine momentum (comparing recent to long-term)
    var momentum = '';
    var momentumColor = '888888';
    if (Math.abs(recentSlopeVal - slopeVal) < 0.002) {
      momentum = '→ Consistent';
      momentumColor = '666666';
    } else if (recentSlopeVal > slopeVal + 0.002) {
      momentum = '↑ Accelerating';
      momentumColor = '228B22';
    } else {
      momentum = '↓ Decelerating';
      momentumColor = 'CC6600';
    }

    inspectorPanel.add(ui.Label('Trend Analysis', { fontWeight: 'bold', fontSize: '11px', margin: '0 0 4px 0' }));

    // 40-year row
    var trend40Panel = ui.Panel({
      widgets: [
        ui.Label('40yr:', { fontSize: '10px', color: '666666' }),
        ui.Label((slopeVal * 1000).toFixed(2) + '×10⁻³/yr', { fontSize: '10px' }),
        ui.Label(slopeClass, { fontSize: '10px', color: slopeColor, fontWeight: 'bold' })
      ],
      layout: ui.Panel.Layout.Flow('horizontal'),
      style: { margin: '0 0 2px 0' }
    });
    inspectorPanel.add(trend40Panel);

    // 10-year row
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

    // Momentum indicator
    inspectorPanel.add(ui.Label(momentum, {
      fontSize: '10px', color: momentumColor, fontStyle: 'italic', margin: '0 0 8px 0'
    }));

    // State transition
    inspectorPanel.add(ui.Label('State Transition', { fontWeight: 'bold', fontSize: '11px', margin: '0 0 4px 0' }));
    inspectorPanel.add(ui.Label(
      '1985: ' + (vegNames[res.start_class] || '?') + '  →  2025: ' + (vegNames[res.end_class] || '?'),
      { fontSize: '11px', margin: '0 0 8px 0' }
    ));

    // Change classification
    inspectorPanel.add(ui.Label('Classification Result', { fontWeight: 'bold', fontSize: '11px', margin: '0 0 4px 0' }));
    if (res.change_class && res.change_class > 0) {
      var changeColors = ['FF00FF', 'FFA500', '7CFC00', '90EE90', '006400', '0000FF', 'FFFF00', '00CED1', 'DDA0DD'];
      var changeColor = changeColors[res.change_class - 1] || 'CCCCCC';
      var lightTextClasses = [3, 4, 7]; // Classes that need dark text
      var changePanel = ui.Panel({
        style: { backgroundColor: '#' + changeColor, padding: '6px 10px', margin: '0 0 4px 0' }
      });
      changePanel.add(ui.Label(classNames[res.change_class] || 'Unknown', {
        fontWeight: 'bold', fontSize: '12px',
        color: lightTextClasses.indexOf(res.change_class) >= 0 ? '000000' : 'FFFFFF'
      }));
      inspectorPanel.add(changePanel);
    } else {
      inspectorPanel.add(ui.Label('No class change detected', { fontSize: '11px', fontStyle: 'italic', color: '888888' }));
    }

    // Canopy Status (unified epoch/projection)
    inspectorPanel.add(ui.Label('Dense Canopy Status', { fontWeight: 'bold', fontSize: '11px', margin: '8px 0 4px 0' }));

    var statusPanel = ui.Panel({ style: { margin: '0 0 4px 0' } });
    var startedDense = res.start_class === 1; // Was Dense at start (1985)
    var nowDense = res.current_ndvi >= 0.6;

    if (res.epoch) {
      // Sparse/Bare → Dense with tracked epoch
      var epochEnd = res.epoch === 2020 ? 2025 : res.epoch + 4;
      statusPanel.add(ui.Label('✓ Established: ' + res.epoch + '-' + epochEnd, {
        fontSize: '11px', color: '228B22', fontWeight: 'bold'
      }));
      statusPanel.add(ui.Label('First reached NDVI ≥0.6 (from sparse/bare)', {
        fontSize: '9px', color: '666666', fontStyle: 'italic'
      }));
    } else if (startedDense && nowDense) {
      // Was Dense, still Dense - established before study
      statusPanel.add(ui.Label('✓ Established: pre-1985', {
        fontSize: '11px', color: '228B22', fontWeight: 'bold'
      }));
      statusPanel.add(ui.Label('Dense canopy throughout study period', {
        fontSize: '9px', color: '666666', fontStyle: 'italic'
      }));
    } else if (!startedDense && nowDense && slopeVal > 0) {
      // Transitioned TO Dense (from Trans) - estimate when it crossed 0.6
      // back-calculate: years ago = (current - 0.6) / slope
      var yearsAgo = (currentNDVI - 0.6) / slopeVal;
      var crossYear = Math.round(2025 - yearsAgo);
      statusPanel.add(ui.Label('✓ Established: ~' + crossYear, {
        fontSize: '11px', color: '228B22', fontWeight: 'bold'
      }));
      statusPanel.add(ui.Label('Estimated crossing from ' + (vegNames[res.start_class] || 'transitional'), {
        fontSize: '9px', color: '666666', fontStyle: 'italic'
      }));
    } else if (!startedDense && nowDense) {
      // Now dense but stable/declining - recently crossed
      statusPanel.add(ui.Label('✓ Established: ~2020s', {
        fontSize: '11px', color: '228B22', fontWeight: 'bold'
      }));
      statusPanel.add(ui.Label('Recently reached threshold', {
        fontSize: '9px', color: '666666', fontStyle: 'italic'
      }));
    } else if (slopeVal > 0) {
      // Gaining - calculate years based on current NDVI and slope
      var yearsNeeded = (0.6 - currentNDVI) / slopeVal;
      if (yearsNeeded <= 0) {
        statusPanel.add(ui.Label('✓ At threshold', {
          fontSize: '11px', color: '228B22', fontWeight: 'bold'
        }));
      } else if (yearsNeeded <= 50) {
        var projYear = 2025 + Math.round(yearsNeeded);
        statusPanel.add(ui.Label('↗ Projected: ~' + projYear, {
          fontSize: '11px', color: '0066CC', fontWeight: 'bold'
        }));
        statusPanel.add(ui.Label('Est. ' + Math.round(yearsNeeded) + ' years to reach NDVI ≥0.6', {
          fontSize: '9px', color: '666666', fontStyle: 'italic'
        }));
      } else {
        statusPanel.add(ui.Label('↗ Projected: ~' + (2025 + Math.round(yearsNeeded)) + ' (slow)', {
          fontSize: '11px', color: '888888'
        }));
        statusPanel.add(ui.Label('~' + Math.round(yearsNeeded) + ' years at current rate', {
          fontSize: '9px', color: '999999', fontStyle: 'italic'
        }));
      }
    } else if (slopeVal < -0.005) {
      // Losing strongly
      statusPanel.add(ui.Label('↘ Declining (no projection)', {
        fontSize: '11px', color: 'CC0000'
      }));
      statusPanel.add(ui.Label('Strong negative trend', {
        fontSize: '9px', color: '999999', fontStyle: 'italic'
      }));
    } else if (slopeVal < 0) {
      // Slight decline
      statusPanel.add(ui.Label('↘ Slight decline', {
        fontSize: '11px', color: 'CC6600'
      }));
      statusPanel.add(ui.Label('Trend negative but below -0.005/yr', {
        fontSize: '9px', color: '999999', fontStyle: 'italic'
      }));
    } else {
      // Truly stable (slope = 0)
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
        title: 'NDVI Trend (1985-2025)',
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
      'Thresholds: Dense≥0.6 | Trans 0.4-0.6 | Sparse 0.2-0.4',
      { fontSize: '9px', color: '888888', margin: '4px 0 0 0' }
    ));
  });
}

Map.onClick(updateInspector);

// 10. EXPORT

Export.image.toDrive({
  image: changeClass.byte(),
  description: 'Export_Change_Classes',
  scale: 30,
  region: roi,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

Export.image.toDrive({
  image: establishmentEpoch.unmask(0).short(),
  description: 'Export_Establishment_Epoch',
  scale: 30,
  region: roi,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});
