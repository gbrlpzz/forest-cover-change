// ============================================================================
// VEGETATION COVER CHANGE DETECTION
// Scientifically Validated Thresholds with Neutral Terminology
// ============================================================================
// This script detects multi-decadal vegetation cover change using a
// regression-first approach with peer-reviewed NDVI thresholds.
//
// SCIENTIFIC SOURCES:
// - NDVI 0.6+ for dense canopy: ResearchGate meta-analyses
// - NDVI 0.2-0.5 sparse/transitional: Copernicus, ISU studies
// - Slope ±0.005/yr threshold: MDPI Kunming study (2000-2020)
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
  max: 6,
  palette: [
    'FF00FF', // 1: Canopy Loss (Magenta)
    'FFA500', // 2: Canopy Thinning (Orange)
    'ADFF2F', // 3: Emerging Biomass (Lime)
    '90EE90', // 4: Canopy Thickening (Light Green)
    '006400', // 5: Canopy Densification (Dark Green)
    '0000FF'  // 6: Canopy Establishment (Blue)
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

// LEGEND
var legend = ui.Panel({
  style: { position: 'bottom-left', padding: '8px 12px', backgroundColor: 'white' }
});

legend.add(ui.Label({
  value: 'Vegetation Change (1985-2025)',
  style: { fontWeight: 'bold', fontSize: '14px', margin: '0 0 8px 0' }
}));

function makeRow(color, name) {
  var colorBox = ui.Label({ style: { backgroundColor: '#' + color, padding: '8px', margin: '0 6px 4px 0' } });
  var desc = ui.Label({ value: name, style: { margin: '0 0 4px 0', fontSize: '11px' } });
  return ui.Panel({ widgets: [colorBox, desc], layout: ui.Panel.Layout.Flow('horizontal') });
}

legend.add(makeRow('FF00FF', 'Canopy Loss'));
legend.add(makeRow('FFA500', 'Canopy Thinning'));
legend.add(makeRow('ADFF2F', 'Emerging Biomass'));
legend.add(makeRow('90EE90', 'Canopy Thickening'));
legend.add(makeRow('006400', 'Canopy Densification'));
legend.add(makeRow('0000FF', 'Canopy Establishment'));

legend.add(ui.Label({ value: 'Thresholds (validated)', style: { fontWeight: 'bold', fontSize: '11px', margin: '8px 0 2px 0' } }));
legend.add(ui.Label({ value: 'Dense: NDVI >= 0.6', style: { fontSize: '10px' } }));
legend.add(ui.Label({ value: 'Transitional: 0.4-0.6', style: { fontSize: '10px' } }));
legend.add(ui.Label({ value: 'Sparse: 0.2-0.4', style: { fontSize: '10px' } }));
legend.add(ui.Label({ value: 'Trend: +/-0.005/yr', style: { fontSize: '10px' } }));

Map.add(legend);

// 9. INSPECTOR

Map.style().set('cursor', 'crosshair');

var classNames = {
  1: 'Canopy Loss',
  2: 'Canopy Thinning',
  3: 'Emerging Biomass',
  4: 'Canopy Thickening',
  5: 'Canopy Densification',
  6: 'Canopy Establishment'
};

var vegNames = { 1: 'Dense Canopy', 2: 'Transitional', 3: 'Sparse', 4: 'Bare' };
var trendNames = { 1: 'Gaining', 2: 'Stable', 3: 'Losing' };

Map.onClick(function (coords) {
  var point = ee.Geometry.Point(coords.lon, coords.lat);

  // Get values at point
  var values = ee.Image.cat([
    changeClass,
    startClass.rename('start_class'),
    endClass.rename('end_class'),
    trendClass.rename('trend_class'),
    slope.rename('slope'),
    intercept.rename('intercept'),
    endNDVI.rename('current_ndvi'),
    establishmentEpoch.rename('epoch'),
    yearsToCanopy.rename('years_proj')
  ]).reduceRegion({
    reducer: ee.Reducer.first(),
    geometry: point,
    scale: 30
  });

  // Create NDVI time series for charting
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

  // Evaluate trend parameters for overlay
  values.evaluate(function (res) {
    // Print summary
    print('═══════════════════════════════════════');
    print('POINT ANALYSIS: ' + coords.lon.toFixed(4) + ', ' + coords.lat.toFixed(4));
    print('═══════════════════════════════════════');

    var slopeVal = res.slope || 0;
    var interceptVal = res.intercept || 0;

    print('Current NDVI: ' + (res.current_ndvi ? res.current_ndvi.toFixed(3) : 'N/A'));
    print('Trend Slope: ' + (slopeVal * 1000).toFixed(3) + ' x10⁻³/yr');
    print('');
    print('1985 State: ' + (vegNames[res.start_class] || 'Unknown'));
    print('2025 State: ' + (vegNames[res.end_class] || 'Unknown'));
    print('40-Year Trend: ' + (trendNames[res.trend_class] || 'Unknown'));

    if (res.change_class && res.change_class > 0) {
      print('');
      print('▶ CHANGE: ' + classNames[res.change_class]);

      if (res.epoch) {
        var epochEnd = res.epoch === 2020 ? 2025 : res.epoch + 4;
        print('  Establishment Epoch: ' + res.epoch + '-' + epochEnd);
      }

      if (res.years_proj && res.years_proj < 50) {
        print('  Projected years to dense canopy: ~' + Math.round(res.years_proj));
      }
    } else {
      print('');
      print('▷ No significant change detected');
    }

    // Create NDVI chart with trend line
    var chart = ui.Chart.feature.byFeature({
      features: ndviSeries,
      xProperty: 'system:time_start',
      yProperties: ['NDVI']
    }).setChartType('ScatterChart')
      .setOptions({
        title: 'NDVI Time Series (Summer Months, 1985-2025)',
        titleTextStyle: { fontSize: 14, bold: true },
        hAxis: {
          title: 'Year',
          format: 'yyyy',
          gridlines: { count: 8 }
        },
        vAxis: {
          title: 'NDVI',
          viewWindow: { min: 0, max: 1 },
          gridlines: { count: 5 }
        },
        pointSize: 3,
        legend: { position: 'bottom' },
        series: {
          0: { color: '333333', pointShape: 'circle' }
        },
        // Add threshold reference lines as annotations
        trendlines: {
          0: {
            type: 'linear',
            color: 'FF0000',
            lineWidth: 2,
            opacity: 0.8,
            showR2: true,
            visibleInLegend: true,
            labelInLegend: 'Linear Trend'
          }
        },
        chartArea: { width: '80%', height: '65%' }
      });

    print(chart);

    // Print threshold reference
    print('───────────────────────────────────────');
    print('Threshold Reference:');
    print('  Dense Canopy: ≥0.6  |  Transitional: 0.4-0.6');
    print('  Sparse: 0.2-0.4     |  Gaining/Losing: ±0.005/yr');
    print('───────────────────────────────────────');
  });
});

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
