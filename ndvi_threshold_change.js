// 1. CONFIGURATION

// Define the period of analysis
var startYear = 1985;
var endYear = 2025;
// NDVI threshold to classify a pixel as forest
var forestThreshold = 0.45;

// Define Region of Interest (ROI). If not explicitly defined, use the map view.
var roi = roi || Map.getBounds(true);

// 2. DATA PROCESSING (Harmonized)
/**
 * Masks clouds/shadows and selects/renames bands for Landsat 5 and 7.
 * Applies scale and offset factors to convert to surface reflectance.
 */
function maskL57(image) {
  var qa = image.select('QA_PIXEL');
  // Cloud and cloud shadow mask (Bits 3 and 4)
  var mask = qa.bitwiseAnd(1 << 3).eq(0).and(qa.bitwiseAnd(1 << 4).eq(0));
  return image.updateMask(mask)
    .select(['SR_B4', 'SR_B3'], ['NIR', 'Red']) // Landsat 5/7 NIR/Red bands
    .multiply(0.0000275).add(-0.2)
    .set('system:time_start', image.get('system:time_start'));
}

/**
 * Masks clouds/shadows and selects/renames bands for Landsat 8 and 9.
 * Applies scale and offset factors to convert to surface reflectance.
 */
function maskL89(image) {
  var qa = image.select('QA_PIXEL');
  // Cloud and cloud shadow mask (Bits 3 and 4)
  var mask = qa.bitwiseAnd(1 << 3).eq(0).and(qa.bitwiseAnd(1 << 4).eq(0));
  return image.updateMask(mask)
    .select(['SR_B5', 'SR_B4'], ['NIR', 'Red']) // Landsat 8/9 NIR/Red bands
    .multiply(0.0000275).add(-0.2)
    .set('system:time_start', image.get('system:time_start'));
}

// Load, filter, and harmonize Landsat collections
var l5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2").filterBounds(roi).map(maskL57);
var l7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2").filterBounds(roi).map(maskL57);
var l8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2").filterBounds(roi).map(maskL89);
var l9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2").filterBounds(roi).map(maskL89);
var fullCollection = l5.merge(l7).merge(l8).merge(l9);

// 3. DEFINE CHANGE CLASSES

// A. START STATE (1985-1989 Average, 5-year window)
// Calculate median NDVI for the baseline period (summer months)
var startNDVI = fullCollection.filterDate('1985-01-01', '1989-12-31')
  .filter(ee.Filter.calendarRange(6, 9, 'month'))
  .median().normalizedDifference(['NIR', 'Red']);
var wasForest = startNDVI.gt(forestThreshold);
var wasOpen = startNDVI.lte(forestThreshold);

// B. END STATE (2021-2025 Average, 5-year window)
// Calculate median NDVI for the final period (summer months)
var endNDVI = fullCollection.filterDate('2021-01-01', '2025-12-31')
  .filter(ee.Filter.calendarRange(6, 9, 'month'))
  .median().normalizedDifference(['NIR', 'Red']);
var isForest = endNDVI.gt(forestThreshold);
var isOpen = endNDVI.lte(forestThreshold);

// C. CHANGE MASKS
// 1. Reforestation (Gain): Was Open -> Is Forest
var reforestationMask = wasOpen.and(isForest);

// 2. Deforestation (Loss): Was Forest -> Is Open
var deforestationMask = wasForest.and(isOpen);

// 4. CALCULATE EPOCH OF REFORESTATION (5-year resolution)

// Define 5-year epochs from 1990 onwards (after baseline period)
var epochs = [
  { start: 1990, end: 1994, label: 1990 },
  { start: 1995, end: 1999, label: 1995 },
  { start: 2000, end: 2004, label: 2000 },
  { start: 2005, end: 2009, label: 2005 },
  { start: 2010, end: 2014, label: 2010 },
  { start: 2015, end: 2019, label: 2015 },
  { start: 2020, end: 2025, label: 2020 }
];

// Create an ImageCollection where each image represents the epoch
// a pixel first crossed the forest threshold.
var epochCollection = ee.ImageCollection.fromImages(
  epochs.map(function (epoch) {
    var startDate = epoch.start + '-01-01';
    var endDate = epoch.end + '-12-31';

    // Get median NDVI for the 5-year epoch (summer months)
    var img = fullCollection
      .filterDate(startDate, endDate)
      .filter(ee.Filter.calendarRange(6, 9, 'month'))
      .median()
      .normalizedDifference(['NIR', 'Red']);

    // Map: 1 if forest, 0 if not. Multiply by epoch label to store the period.
    // .selfMask() removes 0 values (non-forest pixels)
    return img.gt(forestThreshold)
      .multiply(epoch.label)
      .selfMask()
      .toInt()
      .rename('recovery_epoch')
      .set('epoch', epoch.label);
  })
);

// Find the first epoch it crossed the threshold (the minimum epoch value)
var recoveryEpochRaw = epochCollection.min().clip(roi);

// Apply the Reforestation Mask to only show genuine gains (Open -> Forest)
var finalRecoveryMap = recoveryEpochRaw.updateMask(reforestationMask);

// 5. VISUALIZATION

Map.centerObject(roi);

// Layer 1: Deforestation (Loss)
Map.addLayer(deforestationMask.selfMask(),
  { palette: ['FF00FF'] }, // Magenta for Loss
  'Deforestation / Loss of Forest');

// Layer 2: Reforestation (Gain)
// Visualize epoch of recovery using a color gradient (7 distinct epochs)
var vizParams = {
  min: 1990,
  max: 2020,
  palette: [
    '08306b', // 1990-1994: Dark Blue
    '2171b5', // 1995-1999: Blue
    '4eb3d3', // 2000-2004: Light Blue
    '7fcdbb', // 2005-2009: Teal
    'c7e9b4', // 2010-2014: Light Green
    'ffffb2', // 2015-2019: Yellow
    'fd8d3c'  // 2020-2025: Orange
  ]
};
Map.addLayer(finalRecoveryMap, vizParams, 'Reforestation (Epoch of Detection)');


// 6. INSPECTOR

Map.style().set('cursor', 'crosshair');

// Define the action to run on map click
Map.onClick(function (coords) {
  var point = ee.Geometry.Point(coords.lon, coords.lat);

  // Extract values from the final change maps at the clicked point
  var values = ee.Image.cat([
    finalRecoveryMap.rename('recovery_epoch'),
    deforestationMask.rename('is_loss')
  ]).reduceRegion({
    reducer: ee.Reducer.first(),
    geometry: point,
    scale: 30
  });

  // Chart the NDVI time series for the clicked point
  var chart = ui.Chart.image.series({
    imageCollection: fullCollection.select(['NIR', 'Red']),
    region: point,
    reducer: ee.Reducer.median(),
    scale: 30
  }).map(function (img) {
    // Calculate NDVI for each image in the series
    return img.normalizedDifference(['NIR', 'Red']).rename('NDVI');
  }).setOptions({
    title: 'NDVI History (1985-2025)',
    vAxis: { title: 'NDVI', viewWindow: { min: 0, max: 1 } },
    hAxis: { title: 'Year', format: '####' },
    lineWidth: 1,
    pointSize: 2,
    series: { 0: { color: '000000' } }
  });

  // Evaluate the results and print a summary to the console
  values.evaluate(function (res) {
    print('--- Point Analysis ---');

    if (res.recovery_epoch) {
      var epochEnd = res.recovery_epoch === 2020 ? 2025 : res.recovery_epoch + 4;
      print('REFORESTATION DETECTED');
      print('Recovery Epoch: ' + res.recovery_epoch + '-' + epochEnd);
    } else if (res.is_loss === 1) {
      print('DEFORESTATION DETECTED');
      print('Land was Forest in 1985, but is Open today.');
    } else {
      print('STABLE AREA (No change detected)');
    }

    print(chart);
  });
});

// 7. EXPORT CONFIGURATION

// 1. Export the Reforestation Map (Year of Gain)
// Unmask(0) assigns 0 to all "No Reforestation" pixels.
Export.image.toDrive({
  image: finalRecoveryMap.unmask(0).short(), // short() uses 16-bit integer
  description: 'Export_Reforestation_Year',
  scale: 30, // Landsat resolution
  region: roi,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// 2. Export the Deforestation Map (Loss Mask)
// byte() uses 8-bit integer for simple 0/1 data.
Export.image.toDrive({
  image: deforestationMask.unmask(0).byte(),
  description: 'Export_Deforestation_Mask',
  scale: 30,
  region: roi,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});
