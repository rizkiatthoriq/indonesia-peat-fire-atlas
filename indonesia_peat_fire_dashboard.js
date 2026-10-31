/**
 * Indonesia peatland fire activity | Earth Engine Code Editor JavaScript
 * Period: 1 October 2025 through 30 October 2026 UTC.
 * Paste the whole file in a NEW Earth Engine script and press Run.
 * No drawn geometry is needed. Indonesia is the default selection.
 *
 * Descriptive monitoring, not a calibrated fire-risk or smoke-exposure model.
 * The FIRMS source is a raster of ~1 km MODIS detection pixels, not a vector
 * fire-perimeter or an independent-fire count. Confidence >95 is an algorithm
 * score, not a statistical 95% confidence interval.
 *
 * Methods: Giglio et al. (2016), RSE 178:31-41; Giglio et al. (2018), RSE
 * 217:72-85; Funk et al. (2015), Scientific Data 2:150066; Karra et al.
 * (2021), IGARSS; Diggle (1985), JRSS-C 34:138-147.
 */

// --------------------------- Configuration -------------------------------
var CFG = {
  start: '2025-10-01',
  endExclusive: '2026-10-31',
  confidenceGt: 95,
  bandwidthM: 5000,
  kernelRadiusM: 15000,
  highKde: 0.5, // highest display class, detections/km2; descriptive cutoff
  binDays: 14,
  landCoverStart: '2024-01-01',
  landCoverEnd: '2025-01-01',
  colors: ['fff7bc', 'fec44f', 'fe9929', 'cc4c02', '662506'],
  breaks: [0.001, 0.01, 0.05, 0.15, 0.5]
};
var START = ee.Date(CFG.start);
var END = ee.Date(CFG.endExclusive);
var MONTHS = 13;
var LABELS = [];
for (var mi = 0; mi < MONTHS; mi++) {
  LABELS.push(new Date(Date.UTC(2025, 9 + mi, 1))
    .toISOString().slice(0, 7));
}

// ----------------------------- Data --------------------------------------
var provinces = ee.FeatureCollection('FAO/GAUL/2025/level1')
  .filter(ee.Filter.eq('ISO3_CODE', 'IDN'));
// A small bounding box keeps map requests compact. All mapped results are
// additionally restricted by the actual GAUL country mask below.
var indonesiaBox = ee.Geometry.Rectangle([94, -11, 142, 7], null, false);
var countryMask = ee.Image.constant(0).byte()
  .paint(provinces, 1).eq(1);
var peat = ee.Image('projects/sat-io/open-datasets/GLOBAL-PEATLAND-DATABASE')
  .select(0);
var corePeat = peat.eq(1);
var mosaicPeat = peat.eq(2);
var anyPeat = corePeat.or(mosaicPeat);
var fir = ee.ImageCollection('FIRMS')
  .filterBounds(indonesiaBox).filterDate(START, END);
var fireProjection = ee.Image(ee.ImageCollection('FIRMS').first())
  .select('confidence').projection();
var fireZero = ee.Image.constant(0).rename('detections')
  .toInt16().setDefaultProjection(fireProjection);
var rain = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
  .filterBounds(indonesiaBox).filterDate(START, END);
var burned = ee.ImageCollection('MODIS/061/MCD64A1')
  .filterBounds(indonesiaBox).filterDate(START, END);
var landCover = ee.ImageCollection(
  'projects/sat-io/open-datasets/landcover/ESRI_Global-LULC_10m_TS')
  .filterDate(CFG.landCoverStart, CFG.landCoverEnd)
  .filterBounds(indonesiaBox).mosaic().select(0).rename('cover');
landCover = landCover.updateMask(countryMask)
  .updateMask(landCover.neq(10)).updateMask(landCover.neq(0));
var populationSource = ee.ImageCollection(
  'projects/sat-io/open-datasets/WORLDPOP/pop')
  .filter(ee.Filter.stringContains('system:index', 'IDN_POP_2025_'));
var population2025 = ee.Image(populationSource.first())
  .select('population');
// Add population counts, not averages, to the ~1 km fire grid for display
// and proximity summaries. Native 100 m cells are too small to read nationally.
var populationGrid = population2025.updateMask(countryMask)
  .reduceResolution({reducer: ee.Reducer.sum(), maxPixels: 1024})
  .reproject({crs: fireProjection}).rename('people');
var burnedAny = burned.map(burnMask).max().setDefaultProjection(
  ee.Image(burned.first()).select('BurnDate').projection());
var rainfallTotal = rain.select('precipitation').sum();

function makeKernel() {
  return ee.Kernel.gaussian({
    radius: CFG.kernelRadiusM, sigma: CFG.bandwidthM,
    units: 'meters', normalize: true
  });
}
var kernel = makeKernel();

function peatMask(choice) {
  var selected = choice === 'Peat dominated' ? corePeat :
    choice === 'Peat mosaic' ? mosaicPeat : anyPeat;
  return selected.and(countryMask);
}

function detection(image) {
  return image.select('confidence').gt(CFG.confidenceGt)
    .unmask(0).rename('detections').toInt16();
}

// Counting at the native FIRMS grid avoids claiming these are exact point
// locations. Source pixels are masked to peat before KDE, which is then
// restricted to peat again for the map and area summaries.
function detectionCount(start, end, mask) {
  return ee.ImageCollection([fireZero])
    .merge(fir.filterDate(start, end).map(detection))
    .sum().rename('detections')
    .updateMask(mask).unmask(0)
    .reproject({crs: fireProjection});
}

function kde(start, end, mask) {
  return detectionCount(start, end, mask)
    .convolve(kernel)
    .divide(ee.Image.pixelArea().divide(1e6))
    .rename('detections_per_km2')
    .updateMask(mask);
}

function coloredDensity(image) {
  var classes = ee.Image.constant(0).toByte();
  for (var i = 0; i < CFG.breaks.length; i++) {
    classes = classes.where(image.gte(CFG.breaks[i]), i + 1);
  }
  return classes.updateMask(image.gte(CFG.breaks[0]));
}

function burnMask(image) {
  var qa = image.select('QA');
  var valid = qa.bitwiseAnd(1).neq(0)
    .and(qa.rightShift(1).bitwiseAnd(1).neq(0));
  return image.select('BurnDate').gt(0).and(valid)
    .rename('burned').selfMask();
}

function sumIn(image, band, area, scale, crs) {
  var args = {reducer: ee.Reducer.sum(), geometry: area,
    scale: scale, maxPixels: 1e10, tileScale: 16};
  if (crs) args.crs = crs;
  return ee.Dictionary(image.reduceRegion(args)).get(band);
}

function hectares(image, area, scale, crs) {
  return sumIn(ee.Image.pixelArea().divide(10000)
    .rename('ha').updateMask(image).unmask(0),
    'ha', area, scale, crs);
}

// ------------------------- Layout and map --------------------------------
ui.root.clear();
var map = ui.Map();
map.setOptions('ROADMAP');
map.setCenter(117, -2.6, 4);
map.style().set('stretch', 'both');
map.setControlVisibility(true, false, true, true, true, true, false);
ui.root.add(map);
var side = ui.Panel({style: {position: 'top-right', width: '400px',
  height: '60%', padding: '12px', backgroundColor: '#f8faf9'}});
map.add(side);

var ink = '#17332e';
var muted = '#58706a';
function heading(text, size) {
  return ui.Label(text, {fontSize: size || '15px', fontWeight: 'bold',
    color: ink, margin: '12px 0 5px 0'});
}
function note(text) {
  return ui.Label(text, {fontSize: '11px', color: muted,
    margin: '3px 0 7px 0'});
}
function card(title) {
  var value = ui.Label('—', {fontSize: '20px', fontWeight: 'bold',
    color: ink, margin: '0 0 2px 0'});
  var box = ui.Panel({style: {border: '1px solid #dbe5df',
    backgroundColor: 'ffffff', padding: '8px', margin: '4px 0'}});
  box.add(note(title));
  box.add(value);
  side.add(box);
  return value;
}

side.add(ui.Label('INDONESIA / PEATLAND FIRE', {
  fontSize: '11px', fontWeight: 'bold', color: '#35645a',
  margin: '0 0 5px 0'}));
side.add(ui.Label('Where has fire activity concentrated?', {
  fontSize: '22px', fontWeight: 'bold', color: ink,
  margin: '0 0 7px 0'}));
side.add(note('1 Oct 2025–30 Oct 2026 UTC. ' +
  'The map shows smoothed MODIS active-fire detections on mapped peatland.'));

side.add(heading('Explore'));
side.add(note('Province'));
var provinceSelect = ui.Select({items: ['Indonesia'], value: 'Indonesia',
  style: {stretch: 'horizontal'}});
side.add(provinceSelect);
side.add(ui.Button({label: 'Fit selected area on map',
  onClick: function() {
    var selectedName = provinceSelect.getValue() || 'Indonesia';
    if (selectedName === 'Indonesia') map.setCenter(117, -2.6, 4);
    else map.centerObject(selectedGeometry(), 7);
  }, style: {stretch: 'horizontal'}}));
side.add(note('Peatland class'));
var peatSelect = ui.Select({
  items: ['Both mapped classes', 'Peat dominated', 'Peat mosaic'],
  value: 'Both mapped classes', style: {stretch: 'horizontal'}
});
side.add(peatSelect);
side.add(note('Smoothing distance · compare these to test sensitivity'));
var bandwidthSelect = ui.Select({items: ['2 km', '5 km', '10 km'],
  value: '5 km', style: {stretch: 'horizontal'}});
side.add(bandwidthSelect);
side.add(note('Map view'));
var modeSelect = ui.Select({
  items: ['Cumulative since Oct 2025', 'Single month'],
  value: 'Cumulative since Oct 2025',
  style: {stretch: 'horizontal'}
});
side.add(modeSelect);
var monthText = ui.Label('Through 30 Oct 2026',
  {fontWeight: 'bold', color: ink, margin: '7px 0 0 0'});
side.add(monthText);
var monthSlider = ui.Slider({min: 0, max: 12, value: 12, step: 1,
  style: {stretch: 'horizontal'}});
side.add(monthSlider);

side.add(heading('Map key'));
side.add(note('Use these checkboxes to show layers. Ordered quantities ' +
  'use light-to-dark scales; land cover and peat are categories. ' +
  'Rain, land cover and population display one at a time.'));
var legendPanel = ui.Panel();
side.add(legendPanel);

side.add(heading('Selected area'));
side.add(note('Indonesia is the default and includes all mapped peat. ' +
  'Scorecards summarize the full study period; the month ' +
  'slider changes only the map. Burned area ends at its latest source month.'));
var status = note('Calculating Indonesia peatland statistics…');
side.add(status);
var peatAreaCard = card('Mapped peatland area (ha; dominated / mosaic)');
var detectionCard = card('Detection-pixel observations, >95 score');
var highAreaCard = card('Peatland inside highest cumulative KDE class (ha)');
var burnAreaCard = card('Mapped burned peatland (ha; through latest available month)');
var highBurnCard = card('Mapped burned peatland within high KDE class (ha)');
side.add(heading('Residents near peat fire detections'));
var people0to2Card = card('Within 0–2 km (people)');
var people2to5Card = card('2–5 km (people)');
var people5to10Card = card('5–10 km (people)');
var people0to10Card = card('Total within 10 km (people)');
side.add(note('High-activity area is a KDE overlap, not burned area. ' +
  'Distance bands count residents near detected peat fire pixels, ' +
  'including homes outside peatland. These are proximity estimates, ' +
  'not people harmed or exposed to smoke.'));

var chartPanel = ui.Panel({style: {position: 'bottom-right', width: '400px',
  height: '25%', padding: '12px', backgroundColor: '#f8faf9'}});
chartPanel.add(heading('Two-week history'));
chartPanel.add(note('Bars: fire detections. Line: 14-day mean peat ' +
  'rainfall (mm). Missing rainfall remains blank.'));
var fireChartSlot = ui.Panel();
chartPanel.add(fireChartSlot);
map.add(chartPanel);

side.add(heading('Land cover in high-activity class'));
side.add(note('Pre-period Esri annual land cover, 2024. ' +
  'Estimated hectares, sampled at 1 km nationally or 100 m by province. ' +
  'National class areas are coarse.'));
var coverRows = ui.Panel();
side.add(coverRows);
side.add(heading('Land cover in mapped burned area'));
side.add(note('A separate MCD64A1 burned-area estimate, classified ' +
  'using the same pre-period land cover.'));
var burnCoverRows = ui.Panel();
side.add(burnCoverRows);

side.add(heading('Sources and reading'));
side.add(note('Peat: Global Peatland Map 2.0, 1 km (CC BY-NC-SA 4.0). ' +
  'Provinces: FAO GAUL 2025. Fires: NASA FIRMS MODIS NRT, ~1 km. ' +
  'Rain: CHIRPS daily, ~5.6 km. Burned area: MCD64A1, 500 m. ' +
  'Land cover: Esri 2024, 10 m source. Population: WorldPop 2025, 100 m.'));
var freshness = note('Checking source dates…');
side.add(freshness);

var blank = ee.Image(0).selfMask();
var fireLayer = ui.Map.Layer(blank, {}, 'KDE');
var peatLayer = ui.Map.Layer(blank,
  {min: 1, max: 2, palette: ['b8d6c8', 'e0e7dc']},
  'Mapped peatland', true, 0.65);
var rainLayer = ui.Map.Layer(blank,
  {min: 0, max: 3500, palette: ['fff7ec', '9ecae1', '084594']},
  'Rainfall total (mm; available dates)', false, 0.8);
var coverLayer = ui.Map.Layer(blank,
  {min: 1, max: 8, palette: ['1a5bab', '358221', '87d19e',
    'ffdb5c', 'ed022a', 'c4b6aa', 'f2faff', 'c6ad8d']},
  'Land cover (2024)', false, 0.85);
var popLayer = ui.Map.Layer(blank,
  {min: 0, max: 500, palette: ['fff7fb', 'd4b9da', 'c994c7', '980043']},
  'Population (2025; people per ~1 km cell)', false, 0.9);
var burnLayer = ui.Map.Layer(blank, {palette: ['b33422']},
  'Mapped burned area (available months)', false, 0.85);
var detectionLayer = ui.Map.Layer(blank,
  {min: 1, max: 20, palette: ['fee8c8', 'fdbb84', 'e34a33']},
  'Fire detection pixels (>95; full period)', false, 0.9);
var boundaryLayer = ui.Map.Layer(provinces.style({color: '6d3f42',
  fillColor: '00000000', width: 0.5}), {},
  'Indonesia province boundaries', true, 0.85);
map.layers().add(peatLayer);
map.layers().add(rainLayer);
map.layers().add(coverLayer);
map.layers().add(popLayer);
map.layers().add(burnLayer);
map.layers().add(detectionLayer);
map.layers().add(fireLayer);
map.layers().add(boundaryLayer);

var contextBoxes = {};
function layerToggle(label, layer, shown, contextName) {
  var box = ui.Checkbox({label: label, value: shown,
    onChange: function(value) {
      if (contextName) {
        if (value) setContextImage(contextName, selectedCollection());
        else layer.setEeObject(blank);
      }
      layer.setShown(value);
      if (value && contextName) {
        Object.keys(contextBoxes).forEach(function(other) {
          if (other !== contextName) contextBoxes[other].setValue(false);
        });
      }
    },
    style: {margin: '6px 0 1px 0'}});
  legendPanel.add(box);
  return box;
}
function colorMix(a, b, fraction) {
  function channel(i) {
    var x = parseInt(a.substr(i, 2), 16);
    var y = parseInt(b.substr(i, 2), 16);
    var v = Math.round(x + (y - x) * fraction).toString(16);
    return ('0' + v).slice(-2);
  }
  return channel(0) + channel(2) + channel(4);
}
function gradientBar(palette, left, right) {
  var bar = ui.Panel({layout: ui.Panel.Layout.flow('horizontal'),
    style: {margin: '0 0 0 23px'}});
  for (var i = 0; i < 24; i++) {
    var t = i / 23 * (palette.length - 1);
    var j = Math.min(Math.floor(t), palette.length - 2);
    bar.add(ui.Label('', {backgroundColor: '#' + colorMix(
      palette[j], palette[j + 1], t - j),
      width: '10px', height: '11px', margin: '0', padding: '0'}));
  }
  legendPanel.add(bar);
  legendPanel.add(note('     ' + left + '                              ' + right));
}
function categoryKey(color, label) {
  var row = ui.Panel({layout: ui.Panel.Layout.flow('horizontal'),
    style: {margin: '0 0 0 24px'}});
  row.add(ui.Label('', {backgroundColor: '#' + color,
    width: '12px', height: '11px', margin: '3px 5px 0 0'}));
  row.add(note(label));
  legendPanel.add(row);
}
layerToggle('Cumulative / monthly KDE', fireLayer, true);
gradientBar(CFG.colors, '<0.001', '≥0.5 obs/km²');
legendPanel.add(note('     Five fixed classes; darkest is the high-activity class.'));
layerToggle('Mapped peatland', peatLayer, true);
categoryKey('b8d6c8', 'Peat dominated');
categoryKey('e0e7dc', 'Peat in soil mosaic');
contextBoxes.rain = layerToggle('Rainfall total', rainLayer, false, 'rain');
gradientBar(['fff7ec', '9ecae1', '084594'], '0', '≥3,500 mm');
contextBoxes.cover = layerToggle('2024 land cover', coverLayer, false, 'cover');
categoryKey('1a5bab', 'Water');
categoryKey('358221', 'Trees / plantations');
categoryKey('87d19e', 'Flooded vegetation');
categoryKey('ffdb5c', 'Crops');
categoryKey('ed022a', 'Built');
categoryKey('c4b6aa', 'Bare ground');
categoryKey('c6ad8d', 'Rangeland');
contextBoxes.population = layerToggle('Population 2025', popLayer,
  false, 'population');
gradientBar(['fff7fb', 'd4b9da', 'c994c7', '980043'],
  '0', '≥500 people / ~1 km cell');
layerToggle('Mapped burned area', burnLayer, false);
categoryKey('b33422', 'Burned pixel');
layerToggle('Fire detection pixels', detectionLayer, false);
gradientBar(['fee8c8', 'fdbb84', 'e34a33'],
  '1', '≥20 observations');
layerToggle('Indonesia province boundaries', boundaryLayer, true);

// The GIF is a snapshot of the current map extent when the button is pressed.
// A rendered terrain base is used because the interactive Google basemap is
// not an Earth Engine image and cannot be included in a video thumbnail.
var gifPanel = ui.Panel({style: {position: 'bottom-left', width: '320px', height: '320px',
  padding: '7px', backgroundColor: '#f8faf9'}});
gifPanel.add(ui.Label('CUMULATIVE FIRE ACTIVITY / ANIMATION', {
  fontSize: '11px', fontWeight: 'bold', color: ink,
  margin: '0 0 5px 0'}));
gifPanel.add(note('Pan or zoom, then update the GIF for this view.'));
var gifSlot = ui.Panel();
var gifStatus = note('Ready to animate the current map view.');
var gifButton = ui.Button({label: 'Update GIF from current view',
  style: {stretch: 'horizontal'}, onClick: makeGif});
gifPanel.add(gifButton);gifPanel.add(gifStatus);
gifPanel.add(gifSlot);
gifPanel.add(note('Terrain base, province lines, 13 month labels. ' +
  'Ends 30 Oct 2026. Right-click to save.'));
map.add(gifPanel);
var animatedText = require('users/gena/packages:text');
var MONTH_NAMES = ['Oct 2025', 'Nov 2025', 'Dec 2025',
  'Jan 2026', 'Feb 2026', 'Mar 2026', 'Apr 2026',
  'May 2026', 'Jun 2026', 'Jul 2026', 'Aug 2026',
  'Sep 2026', 'Oct 2026'];

function makeGif() {
  var bounds = map.getBounds(); // west, south, east, north in degrees
  var view = ee.Geometry.Rectangle(bounds, null, false);
  var selected = selectedCollection();
  var mask = peatMask(peatSelect.getValue());
  var centreLat = (bounds[1] + bounds[3]) / 2;
  var metresPerPixel = (bounds[2] - bounds[0]) * 111320 *
    Math.cos(centreLat * Math.PI / 180) / 250;
  var textScale = Math.max(100, metresPerPixel * 2);
  var base = ee.Image('NOAA/NGDC/ETOPO1').select('bedrock')
    .visualize({min: -100, max: 800,
      palette: ['d7e9ed', 'e7eee6', 'd7e4d5', 'd1d5bd', 'e7e0c9']});
  var outlines = provinces.filterBounds(view).style({
    color: '516a62', fillColor: '00000000', width: 1});
  var frames = [];
  for (var i = 0; i < MONTHS; i++) {
    var end = i === 12 ? END : START.advance(i + 1, 'month');
    var activity = coloredDensity(kde(START, end, mask))
      .clipToCollection(selected)
      .visualize({min: 1, max: 5, palette: CFG.colors, opacity: 0.85});
    var frame = base.blend(activity).blend(outlines)
      .set('label', MONTH_NAMES[i]);
    frames.push(animatedText.annotateImage(frame,
      {forceRgbOutput: true}, view, [{
        position: 'left', offset: '5%', margin: '4%',
        property: 'label', scale: textScale,
        fontType: 'Consolas', fontSize: 18,
        textColor: '17332e', outlineColor: 'ffffff',
        outlineWidth: 3
      }]));
  }
  gifSlot.clear();
  gifStatus.setValue('Loading 13 frames; allow a few minutes.');
  gifSlot.add(ui.Thumbnail({image: ee.ImageCollection.fromImages(frames),
    params: {region: view, dimensions: '250x140',
      framesPerSecond: 1.4, crs: 'EPSG:3857', format: 'gif'},
    style: {width: '250px', height: '140px'}}));
}

// ---------------------- On-demand computation ----------------------------
var selectionSerial = 0;
function selectedGeometry() {
  var name = provinceSelect.getValue() || 'Indonesia';
  return name === 'Indonesia' ? indonesiaBox :
    provinces.filter(ee.Filter.eq('GAUL1_NAME', name)).geometry();
}
function selectedCollection() {
  var name = provinceSelect.getValue() || 'Indonesia';
  return name === 'Indonesia' ? provinces :
    provinces.filter(ee.Filter.eq('GAUL1_NAME', name));
}
function setContextImage(name, region) {
  if (name === 'rain') rainLayer.setEeObject(
    rainfallTotal.updateMask(countryMask).clipToCollection(region));
  if (name === 'cover') {
    var remappedCover = landCover.remap([1, 2, 4, 5, 7, 8, 9, 11],
      [1, 2, 3, 4, 5, 6, 7, 8]);
    coverLayer.setEeObject(remappedCover.clipToCollection(region));
  }
  if (name === 'population') popLayer.setEeObject(
    populationGrid.clipToCollection(region));
}
function displayedDates() {
  var n = Math.round(monthSlider.getValue());
  var start = START.advance(n, 'month');
  var end = n === 12 ? END : start.advance(1, 'month');
  return {start: modeSelect.getValue() === 'Single month' ? start : START,
    end: end, index: n};
}
function refreshMap() {
  var dates = displayedDates();
  var mask = peatMask(peatSelect.getValue());
  var region = selectedCollection();
  var density = kde(dates.start, dates.end, mask);
  var name = modeSelect.getValue() === 'Single month' ?
    'Activity in ' : 'Accumulated to ';
  fireLayer.setEeObject(coloredDensity(density).clipToCollection(region));
  fireLayer.setVisParams({min: 1, max: 5, palette: CFG.colors});
  fireLayer.setName(name + LABELS[dates.index]);
  peatLayer.setEeObject(peat.updateMask(mask).clipToCollection(region));
  // Context maps cover the administrative area and are requested only when
  // visible. Hidden national rasters otherwise delay the fire map on startup.
  Object.keys(contextBoxes).forEach(function(contextName) {
    if (contextBoxes[contextName].getValue())
      setContextImage(contextName, region);
  });
  burnLayer.setEeObject(burnedAny.updateMask(mask).selfMask()
    .clipToCollection(region));
  var fullCount = detectionCount(START, END, mask);
  detectionLayer.setEeObject(fullCount.updateMask(fullCount.gt(0))
    .clipToCollection(region));
  boundaryLayer.setEeObject(region.style({color: '6d3f42',
    fillColor: '00000000', width: 0.5}));
  monthText.setValue((modeSelect.getValue() === 'Single month' ?
    'During ' : 'Through ') + LABELS[dates.index] +
    (dates.index === 12 ? ' · through 30 Oct' : ''));
}

function formatNumber(value, unit) {
  if (value === null || value === undefined) return 'Unavailable';
  return Math.round(value).toLocaleString('en-US') + (unit || '');
}
function setCard(label, value, unit, serial) {
  if (serial === selectionSerial) label.setValue(formatNumber(value, unit));
}
function cardResult(label, value, error, unit, serial, description) {
  if (serial !== selectionSerial) return;
  if (error) {
    label.setValue('Error — see Console');
    print(description + ' failed', error);
  } else {
    setCard(label, value, unit, serial);
  }
}
function updateCards(area, mask, serial) {
  var coreArea = hectares(corePeat.and(countryMask), area, 1000, null);
  var mosaicArea = hectares(mosaicPeat.and(countryMask), area, 1000, null);
  ee.Dictionary({core: coreArea, mosaic: mosaicArea}).evaluate(function(v, e) {
    if (serial !== selectionSerial) return;
    if (e) {
      peatAreaCard.setValue('Error — see Console');
      print('Mapped peat area failed', e);
    } else peatAreaCard.setValue(formatNumber(v.core) + ' / ' +
      formatNumber(v.mosaic));
  });

  var count = detectionCount(START, END, mask).updateMask(mask);
  sumIn(count, 'detections', area, 1000, fireProjection)
    .evaluate(function(v, e) {
      cardResult(detectionCard, v, e, '', serial, 'FIRMS count');
    });

  var high = kde(START, END, mask).gte(CFG.highKde).selfMask();
  hectares(high, area, 1000, fireProjection).evaluate(function(v, e) {
    cardResult(highAreaCard, v, e, ' ha', serial, 'High KDE area');
  });

  var burnAll = burnedAny.updateMask(mask);
  hectares(burnAll, area, 500, null).evaluate(function(v, e) {
    cardResult(burnAreaCard, v, e, ' ha', serial, 'Burned peat area');
  });
  var burnedFraction = burnAll.unmask(0)
    .reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 64})
    .reproject({crs: fireProjection});
  sumIn(ee.Image.pixelArea().divide(10000).multiply(burnedFraction)
    .rename('ha').updateMask(high), 'ha', area, 1000, fireProjection)
    .evaluate(function(v, e) {
      cardResult(highBurnCard, v, e, ' ha', serial,
        'Burned area within high KDE');
    });

  // Sum WorldPop people per cell into the FIRMS grid before masking by
  // distance. The source detections are on peat, but residents may live
  // on surrounding non-peat land. Bands are exclusive and sum to 0–10 km.
  var peopleOnFireGrid = populationGrid;
  var sourcePixels = detectionCount(START, END, mask).gt(0)
    .clipToCollection(selectedCollection()).unmask(0);
  var distanceM = sourcePixels.distance(
    ee.Kernel.euclidean({radius: 10000, units: 'meters'}));
  var near2 = distanceM.lte(2000);
  var near5 = distanceM.gt(2000).and(distanceM.lte(5000));
  var near10 = distanceM.gt(5000).and(distanceM.lte(10000));
  var populationBands = ee.Image.cat([
    peopleOnFireGrid.updateMask(near2).rename('people_0_2'),
    peopleOnFireGrid.updateMask(near5).rename('people_2_5'),
    peopleOnFireGrid.updateMask(near10).rename('people_5_10')
  ]);
  populationBands.reduceRegion({reducer: ee.Reducer.sum(),
    geometry: area, scale: 1000, crs: fireProjection,
    maxPixels: 1e10, tileScale: 16}).evaluate(function(v, e) {
      if (serial !== selectionSerial) return;
      if (e || !v) {
        [people0to2Card, people2to5Card, people5to10Card,
          people0to10Card].forEach(function(label) {
            label.setValue('Error — see Console');
          });
        print('WorldPop distance bands failed', e);
        return;
      }
      var p2 = v.people_0_2 || 0;
      var p5 = v.people_2_5 || 0;
      var p10 = v.people_5_10 || 0;
      people0to2Card.setValue(formatNumber(p2));
      people2to5Card.setValue(formatNumber(p5));
      people5to10Card.setValue(formatNumber(p10));
      people0to10Card.setValue(formatNumber(p2 + p5 + p10));
    });

  var cover = landCover;
  var coverArea = ee.Image.pixelArea().divide(10000).rename('ha')
    .addBands(cover).updateMask(high);
  function showCover(image, slot) {
    var reportingScale = provinceSelect.getValue() === 'Indonesia' ||
      !provinceSelect.getValue() ? 1000 : 100;
    image.reduceRegion({reducer: ee.Reducer.sum().group({
      groupField: 1, groupName: 'class'}), geometry: area,
      scale: reportingScale,
      maxPixels: 1e10, tileScale: 16}).evaluate(function(v, e) {
      if (serial !== selectionSerial) return;
      slot.clear();
      if (e || !v || !v.groups || !v.groups.length) {
        if (e) print('Land-cover hectares failed', e);
        slot.add(note('No classified pixels, or summary unavailable.'));
        return;
      }
      var names = {1: 'Water', 2: 'Trees / plantations',
        4: 'Flooded vegetation', 5: 'Crops', 7: 'Built',
        8: 'Bare ground', 9: 'Snow / ice', 11: 'Rangeland'};
      v.groups.sort(function(a, b) { return b.sum - a.sum; });
      v.groups.forEach(function(g) {
        slot.add(ui.Label((names[g.class] || 'Other / unclassified') + '  ·  ' +
          formatNumber(g.sum, ' ha'), {fontSize: '12px', color: ink,
          margin: '2px 0'}));
      });
    });
  }
  showCover(coverArea, coverRows);
  showCover(ee.Image.pixelArea().divide(10000).rename('ha')
    .addBands(cover).updateMask(burnAll), burnCoverRows);
}

function updateCharts(area, mask, serial) {
  // Evaluate fire and rainfall separately within each 14-day bin to avoid
  // Earth Engine's concurrent-aggregation limit on a national chart.
  var startMs = Date.parse(CFG.start + 'T00:00:00Z');
  var endMs = Date.parse(CFG.endExclusive + 'T00:00:00Z');
  var dayMs = 86400000;
  var n = Math.ceil((endMs - startMs) / (CFG.binDays * dayMs));
  var rows = [];
  var missingRainBins = 0;
  fireChartSlot.clear();
  var chartProgress = note('Calculating 14-day chart: 0 / ' + n);
  fireChartSlot.add(chartProgress);

  function next(i) {
    if (serial !== selectionSerial) return;
    if (i === n) {
      var series = ee.FeatureCollection(rows);
      fireChartSlot.clear();
      if (missingRainBins) fireChartSlot.add(note(
        'Rainfall could not be calculated for ' + missingRainBins +
        ' bin(s); those line segments are blank. See Console.'));
      fireChartSlot.add(ui.Chart.feature.byFeature(series, 'label',
        ['detections', 'rainfall_mm'])
        .setChartType('ComboChart').setOptions({
          title: 'Fire detections and rainfall · 14-day bins',
          legend: {position: 'bottom'},
          chartArea: {left: 55, right: 55, top: 24, bottom: 55},
          hAxis: {slantedText: true, slantedTextAngle: 60,
            showTextEvery: 4, textStyle: {fontSize: 8}},
          height: 175,
          colors: ['#a3441b', '#427c94'],
          series: {0: {type: 'bars', targetAxisIndex: 0},
            1: {type: 'line', targetAxisIndex: 1, lineWidth: 2}},
          vAxes: {0: {title: 'Detections'}, 1: {title: 'Rainfall (mm)'}}
        }));
      return;
    }
    var fromMs = startMs + i * CFG.binDays * dayMs;
    var tillMs = Math.min(fromMs + CFG.binDays * dayMs, endMs);
    var from = ee.Date(new Date(fromMs).toISOString());
    var till = ee.Date(new Date(tillMs).toISOString());
    var days = (tillMs - fromMs) / dayMs;
    var rainPart = rain.filterDate(from, till).select('precipitation');
    var nRainImages = rainPart.size();
    var rainSum = ee.Image(ee.Algorithms.If(nRainImages.gt(0),
      rainPart.sum(),
      ee.Image.constant(0).toFloat().rename('precipitation').selfMask()))
      .updateMask(mask);
    var mm = ee.Dictionary(rainSum.reduceRegion({
      reducer: ee.Reducer.mean(), geometry: area, scale: 5566,
      maxPixels: 1e8, tileScale: 16})).get('precipitation');
    var hits = sumIn(detectionCount(from, till, mask)
      .updateMask(mask), 'detections', area, 1000, fireProjection);
    function finish(hitsValue, rainValue) {
      rows.push(ee.Feature(null, {
        label: new Date(fromMs).toLocaleDateString('en-GB',
          {day: '2-digit', month: 'short', year: '2-digit',
            timeZone: 'UTC'}),
        detections: hitsValue,
        rainfall_mm: rainValue
      }));
      chartProgress.setValue('Calculating 14-day chart: ' +
        (i + 1) + ' / ' + n);
      next(i + 1);
    }
    ee.Dictionary({detections: ee.Algorithms.If(hits, hits, 0),
      rainDays: nRainImages}).evaluate(function(v, error) {
        if (serial !== selectionSerial) return;
        if (error) {
          chartProgress.setValue('Chart could not be calculated. See Console.');
          print('14-day fire bin ' + (i + 1) + ' failed', error);
          return;
        }
        if (v.rainDays < days) {
          finish(v.detections, null);
        } else {
          mm.evaluate(function(mmValue, rainError) {
            if (serial !== selectionSerial) return;
            if (rainError) {
              print('14-day rainfall bin ' + (i + 1) + ' failed', rainError);
              missingRainBins++;
              finish(v.detections, null);
              return;
            }
            finish(v.detections, mmValue);
          });
        }
      });
  }
  next(0);
}

function refreshDetails(recomputeCharts) {
  selectionSerial++;
  var serial = selectionSerial;
  var name = provinceSelect.getValue() || 'Indonesia';
  status.setValue('Calculating ' + name + '…');
  coverRows.clear();
  burnCoverRows.clear();
  [peatAreaCard, detectionCard, highAreaCard, burnAreaCard, highBurnCard,
    people0to2Card, people2to5Card, people5to10Card,
    people0to10Card].forEach(function(label) {
      label.setValue('Calculating…');
    });
  var area = selectedGeometry();
  var mask = peatMask(peatSelect.getValue());
  if (name === 'Indonesia') map.setCenter(117, -2.6, 4);
  else map.centerObject(area, 7);
  ee.Dictionary(mask.rename('peat').reduceRegion({
    reducer: ee.Reducer.sum(), geometry: area, scale: 1000,
    maxPixels: 1e8, tileScale: 4})).get('peat', 0)
    .evaluate(function(pixels, error) {
      if (serial !== selectionSerial) return;
      if (error) {
        status.setValue('Peat coverage could not be checked for ' + name + '.');
        print('Peat coverage check failed', error);
        return;
      }
      if (!pixels) {
        status.setValue('No mapped pixels of the selected peat class ' +
          'in ' + name + '. Try the other peat class.');
        [peatAreaCard, detectionCard, highAreaCard, burnAreaCard,
          highBurnCard, people0to2Card, people2to5Card,
          people5to10Card, people0to10Card].forEach(function(label) {
            label.setValue('—');
          });
        fireChartSlot.clear();
        coverRows.clear();
        burnCoverRows.clear();
        return;
      }
      status.setValue(name + ' · peat class: ' + peatSelect.getValue() +
        '. Area card separates dominated and mosaic pixels.');
      updateCards(area, mask, serial);
      if (recomputeCharts !== false) updateCharts(area, mask, serial);
    });
}

provinceSelect.onChange(function() {
  refreshMap(); refreshDetails();
});
peatSelect.onChange(function() { refreshMap(); refreshDetails(); });
bandwidthSelect.onChange(function(value) {
  CFG.bandwidthM = parseInt(value, 10) * 1000;
  CFG.kernelRadiusM = 3 * CFG.bandwidthM;
  kernel = makeKernel();
  refreshMap();
  refreshDetails(false);
});
modeSelect.onChange(refreshMap);
monthSlider.onChange(refreshMap);
refreshMap();
refreshDetails();

ee.Image.constant(1).rename('peat').updateMask(anyPeat.and(countryMask))
  .reduceRegions({collection: provinces, reducer: ee.Reducer.count(),
    scale: 1000, tileScale: 8})
  .filter(ee.Filter.gt('count', 0))
  .aggregate_array('GAUL1_NAME').distinct().sort()
  .evaluate(function(names, error) {
    if (error || !names) {
      status.setValue('Peatland province list unavailable; see Console.');
      if (error) print('Peatland province filter failed', error);
      return;
    }
    names.forEach(function(name) {
      if (name !== 'Indonesia') provinceSelect.items().add(name);
    });
  });

ee.Dictionary({fir: fir.aggregate_max('system:time_start'),
  rain: rain.aggregate_max('system:time_start'),
  burn: burned.aggregate_max('system:time_start'),
  populationImages: populationSource.size()})
  .evaluate(function(v, error) {
    if (error || !v) {
      freshness.setValue('Source dates could not be checked.');
      return;
    }
    function date(x) {
      return x ? new Date(x).toISOString().slice(0, 10) : 'none';
    }
    freshness.setValue('Latest in period · FIRMS: ' + date(v.fir) +
      '; CHIRPS: ' + date(v.rain) + '; burned area: ' + date(v.burn) +
      '. WorldPop 2025 images matched: ' + v.populationImages + '.');
  });

print('Dashboard version', 'v4 Indonesia peatland fire activity');
print('Scientific note', 'KDE and its high class indicate relative ' +
  'concentration of repeated detection-pixel observations, not burn severity, ' +
  'fire probability, ecological damage or smoke exposure.');