/**
 * @license
 * Copyright 2020 Justin Braaten
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */



// Ideas for correcting sensors: https://ieeexplore.ieee.org/abstract/document/9093966

var msslib = require('users/jstnbraaten/modules:msslib/msslib.js');
var ltgee = require('users/emaprlab/public:Modules/LandTrendr.js'); 
var animation = require('users/gena/packages:animation');



/**
* Returns a filtered TM WRS-2 T1 surface reflectance image collection.
* 
* @param {ee.Geometry | ee.Feature} aoi An ee.Filter to filter TM image collection.
* 
* @return {ee.ImageCollection} An MSS WRS-2 image collection filtered by
*     bounds and quality.
*/
function getTmWrs2Col(aoi){
  var tm4 = ee.ImageCollection("LANDSAT/LT04/C01/T1_SR")
    .filterBounds(aoi);
  var tm5 = ee.ImageCollection("LANDSAT/LT05/C01/T1_SR")
    .filterBounds(aoi);
  return tm4.merge(tm5);
}
exports.getTmWrs2Col = getTmWrs2Col;




/**
* Add unique path, row, orbit ID as image property for joining TM and MSS collections.
* 
* @param {ee.Image} tmWrs2Col A TM image collection.
* @param {ee.Image} mssWrs2Col A MSS image collection.
* 
* @return {ee.ImageCollection} An image collection ____WAH_____.
*/ 
function coincidentTmMssCol(tmWrs2Col, mssWrs2Col){
  var filter = ee.Filter.equals({leftField: 'imgID', rightField: 'imgID'});
  var join = ee.Join.saveFirst('coincidentTmMss');
  return ee.ImageCollection(join.apply(tmWrs2Col, mssWrs2Col, filter));
}
exports.coincidentTmMssCol = coincidentTmMssCol;




/**
* Add unique path, row, orbit ID as image property for joining TM and MSS collections.
* 
* @param {ee.Image} img A TM or MSS image.
* 
* @return {ee.ImageCollection} A copy of the input image with an 'imgID'
*     property added to the image describing the unique path, row, orbit.
*/ 
function addTmToMssJoinId(img){
  //return col.map(function(img) {
    var date = ee.Image(img).date();
    var year = ee.Algorithms.String(date.get('year'));
    var doy = ee.Algorithms.String(date.getRelative('day', 'year'));
    var path = ee.Algorithms.String(img.getNumber('WRS_PATH').toInt());
    var row = ee.Algorithms.String(img.getNumber('WRS_ROW').toInt());
    var yearDoy = year.cat(doy).cat(path).cat(row);
    return img.set({'imgID': yearDoy,
      'path': path,
      'row': row
    });
  //});
}
exports.addTmToMssJoinId = addTmToMssJoinId;




/**
 * Returns the footprint of an image as a ee.Geometry.Polygon. 
 * 
 * @param {ee.Image} img The image to get the footprint for.
 * 
 * @return {ee.Geometry.Polygon} The ee.Geometry.Polygon representation of
 *     an image's footprint.
 */
function getFootprint(img){
  return ee.Geometry.Polygon(ee.Geometry(img.get('system:footprint')).coordinates())}
exports.getFootprint = getFootprint;



/**
 * Generates an ee.Filter for filtering MSS and TM image collection by
 *     intersection with a given geometry. 
 * 
 * @param {ee.Geometry | ee.Feature} aoi Area of interest to filter collection to.
 *     Include images less than given value.
 * 
 * @return {ee.Filter} A filter to be passed as an argument to the .filter()
 *     ee.ImageCollection method.
 */
function filterBounds(aoi) {
  return ee.Filter.bounds(aoi);
}
exports.filterBounds = filterBounds;


/**
 * Returns a cloud and cloud shadow mask from CFmask. 
 * @param {ee.Image} img Landsat SR image.
 * @return {ee.Image} A 0/1 mask image to be used with .updateMask().
 */
function getCfmask(img) {
  var cloudShadowBitMask = 1 << 3;
  var cloudsBitMask = 1 << 5;
  var qa = img.select('pixel_qa');
  var mask = qa.bitwiseAnd(cloudShadowBitMask)
    .eq(0)
    .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
  return mask;
}
exports.getCfmask = getCfmask;


function applyCfmask(img) {
  var mask = getCfmask(img);
  return img.updateMask(mask);
}
exports.getCfmask = getCfmask;



// #############################################################################
// ### Process steps ###
// #############################################################################


/**
 * Display the series of WRS-1 images for a given WRS-1 granule.
 */
function viewWrs1Col(params) {
  print('Displaying WRS-1 images to the console');
  var granuleGeom = msslib.getWrs1GranuleGeom(params.wrs1);
  params.aoi = ee.Geometry(granuleGeom.get('centroid'));
  params.wrs = '1';
  var mssDnCol = msslib.getCol(params)
    .filter(ee.Filter.eq('pr', params.wrs1));
  msslib.viewThumbnails(mssDnCol);
}
exports.viewWrs1Col = viewWrs1Col;

/**
 * Display the WRS-1 grid to the map.
 */
function wrs1GranuleSelector() {
  var wrs1Granules = ee.FeatureCollection('users/jstnbraaten/wrs/wrs1_descending_land');
  Map.addLayer(wrs1Granules, {color: 'grey'}, null, null, 0.5);
  
  var message = ui.Label({value: 'Click granules to print ID. Wait patiently after clicking. Repeat as needed.',
    style: {position: 'top-center'}});
  var holder = ui.Panel({style: {width: '170px', height: '220px', position: 'top-left'}});
  var label = ui.Label({value: 'WRS-1 Granule ID(s):'});
  var ids = ui.Panel({style: {backgroundColor: '#DCDCDC'}});
  holder.add(label);
  holder.add(ids);
  
  Map.add(message);
  Map.add(holder);
  Map.style().set('cursor', 'crosshair');
  Map.onClick(function(e) {
    ids.clear();
    var nLayers = Map.layers().length();
    for(var i=0; i < nLayers-1; i++) {
      Map.layers().remove(Map.layers().get(1));
    }
    
    var point = ee.Geometry.Point(e.lon, e.lat);
    var joinFilter = ee.Filter.intersects({leftField: '.geo', rightField: '.geo', maxError: 500});
    var join = ee.Join.simple();
    var intersectingFeatures = join.apply(wrs1Granules, ee.FeatureCollection(point), joinFilter);
  
    intersectingFeatures.toList(intersectingFeatures.size()).evaluate(function(fList) {
      var colors = ['red', 'blue', 'green', 'yellow', 'orange', 'pink', 'purple'];
      for(var i in fList) {
        var f = ee.Feature(fList[i]);
        var outline = ee.Image().byte().paint({
          featureCollection: ee.FeatureCollection(f),
          color: 1,
          width: 3
        });
        var pr = f.get('PR').getInfo();
        var title = pr + ' ' + colors[i];
        Map.addLayer(outline, {palette: colors[i]}, title);
        ids.add(ui.Label({value: pr, style: {color: colors[i], backgroundColor: '#DCDCDC'}}));
      }
    });
  });
}
exports.wrs1GranuleSelector = wrs1GranuleSelector;



// #############################################################################
// ### Reference prep ###
// #############################################################################

/**
 * calculate the medoid of a collection. 
 * 
 */
function getMedoid(col, bands) {
  col = col.select(bands);
  var median = col.median();
  var difFromMedian = col.map(function(img) {
    var dif = ee.Image(img).subtract(median).pow(ee.Image.constant(2));
    return dif.reduce(ee.Reducer.sum())
      .addBands(img);
  });
  var bandNames = difFromMedian.first().bandNames();
  var len = bandNames.length();
  var bandsPos = ee.List.sequence(1, len.subtract(1));
  var bandNamesSub = bandNames.slice(1);
  return difFromMedian.reduce(ee.Reducer.min(len)).select(bandsPos, bandNamesSub);
}
exports.getMedoid = getMedoid;

function getRefImg(params) {
  var granuleGeoms = msslib.getWrs1GranuleGeom(params.wrs1);
  var centroid = ee.Geometry(granuleGeoms.get('centroid'));
  var bounds = ee.Geometry(granuleGeoms.get('bounds'));

  var refCol = msslib.getCol({
    aoi: bounds,
    wrs: '2',
    yearRange: [1983, 1987], // Use five early years - want good coverage, but near to MSS WRS-1 window.
    doyRange: params.doyRange,
  }).map(addTmToMssJoinId);

  var tmCol = getTmWrs2Col(bounds)
    .filterDate('1983-01-01', '1988-01-01')
    .map(addTmToMssJoinId);
  
  var coincident = coincidentTmMssCol(refCol, tmCol)
    .map(function(img) {
      var mask = getCfmask(ee.Image(img.get('coincidentTmMss')));
      var imgToa = msslib.calcToa(img);
      return imgToa.updateMask(mask);
    });

  return msslib.addTc(msslib.addNdvi(getMedoid(coincident, ['green', 'red', 'red_edge', 'nir'])))
    .select(['green', 'red', 'nir', 'ndvi', 'tcb', 'tcg', 'tca'])  // TODO: scale these data to get them to int16
    .set('bounds', bounds);
}
exports.getRefImg = getRefImg;


function exportMssRefImg(params) {
  print('Preparing reference image export task, please wait');
  var refImg = getRefImg(params);
  Export.image.toAsset({
    image: refImg,
    description: 'MSS-reference-image',
    assetId: params.baseDir + '/ref',
    region: ee.Geometry(refImg.get('bounds')),
    scale: 60,
    crs: params.crs,
    maxPixels: 1e13
  });  
}
exports.exportMssRefImg = exportMssRefImg;






// #############################################################################




/**
 * Returns an example Tm image. 
 * 
 * @return {ee.Image} Example TM image.
 */
function exampleTmImg() {
  return ee.Image('LANDSAT/LT05/C01/T1_SR/LT05_045029_19840728'); 
}
exports.exampleTmImg = exampleTmImg;






// Example AOIs
var wrs2045029 = ee.Geometry.Point([-121.454, 44.47]);
exports.wrs2045029 = wrs2045029;




// #############################################################################
// ### Process MSS WRS-1 images ###
// #############################################################################


// Create a new image that is the concatenation of three images: a constant,
// the SWIR1 band, and the SWIR2 band.

/**
 * Have found that the best regression is robustLinear on entire image
 * with scale set as 60.
 * Also tried: robustLinear, scale 300
 * linear, scale 300
 * stratified sample based on image segmentation k-means - worst - could be poor sampling
 * Yet to try using all bands to predicit a given band - multiple regression
 * Yest to try using different threshold and scale for masking dif to ref in `correctMssImg` 
 */
function calcRegression(xImg, yImg, xBand, yBand, aoi, scale) {
  var constant = ee.Image(1);
  var xVar = xImg.select(xBand);
  var yVar = yImg.select(yBand);
  var imgRegress = ee.Image.cat(constant, xVar, yVar);
  
  var linearRegression = imgRegress.reduceRegion({
    reducer: ee.Reducer.robustLinearRegression({
      numX: 2,
      numY: 1
    }),
    geometry: aoi,
    scale: scale,
    maxPixels: 1e13
  });

  var coefList = ee.Array(linearRegression.get('coefficients')).toList();
  var intercept = ee.List(coefList.get(0)).get(0);
  var slope = ee.List(coefList.get(1)).get(0);
  var rmse = ee.Array(linearRegression.get('residuals')).toList().get(0);
  return ee.Dictionary({slope: slope, intercept: intercept, rmse: rmse});
}
exports.calcRegression = calcRegression;

// Function to apply correction to reference image.
function applyCoef(img, band, coef) {
  coef = ee.Dictionary(coef);
  return img.select(band)
    .multiply(ee.Image.constant(coef.getNumber('slope')))
    .add(ee.Image.constant(coef.getNumber('intercept')));
}
exports.applyCoef = applyCoef;


function getSampleImg(img, ref, band) {
  var dif = img.select(band)
    .subtract(ref.select(band)).rename('dif');
  
  var difThresh = dif.reduceRegion({
    reducer: ee.Reducer.percentile({
      percentiles: [40, 60],
      maxRaw: 1000000,
      maxBuckets: 1000000,
      minBucketWidth: 0.00000000001
    }),
    geometry: img.geometry(),
    scale: 60,
    maxPixels: 1e13
  });
  
  var mask = dif.gt(difThresh.getNumber('dif_p40'))
    .and(dif.lt(difThresh.getNumber('dif_p60')));

  return img.updateMask(mask);
}
exports.getSampleImg = getSampleImg;


// Function to make normalization function.
function correctMssImg(img) {
  var ref = ee.Image(img.get('ref_img'));
  var granuleGeoms = msslib.getWrs1GranuleGeom(img.getString('pr'));
  var granule = ee.Feature(granuleGeoms.get('granule')).geometry(); 
  
  // // ** Use three bands for mask
  // var allMask = getSampleImg(img, ref, 'green').mask().multiply(
  //   getSampleImg(img, ref, 'red').mask()).multiply(
  //   getSampleImg(img, ref, 'nir').mask());
  // var sampleImg = img.updateMask(allMask);
  
  // var greenCoef = calcRegression(sampleImg, ref, 'green', 'green', granule, 60);
  // var redCoef = calcRegression(sampleImg, ref, 'red', 'red', granule, 60);
  // var nirCoef = calcRegression(sampleImg, ref, 'nir', 'nir', granule, 60);
  // var ndviCoef = calcRegression(sampleImg, ref, 'ndvi', 'ndvi', granule, 60);
  // var tcbCoef = calcRegression(sampleImg, ref, 'tcb', 'tcb', granule, 60);
  // var tcgCoef = calcRegression(sampleImg, ref, 'tcg', 'tcg', granule, 60);
  // var tcaCoef = calcRegression(sampleImg, ref, 'tca', 'tca', granule, 60);
  // // ** Use three bands for mask
  
  var greenCoef = calcRegression(getSampleImg(img, ref, 'green'), ref, 'green', 'green', granule, 60);
  var redCoef = calcRegression(getSampleImg(img, ref, 'red'), ref, 'red', 'red', granule, 60);
  var nirCoef = calcRegression(getSampleImg(img, ref, 'nir'), ref, 'nir', 'nir', granule, 60);
  var ndviCoef = calcRegression(getSampleImg(img, ref, 'ndvi'), ref, 'ndvi', 'ndvi', granule, 60);
  var tcbCoef = calcRegression(getSampleImg(img, ref, 'tcb'), ref, 'tcb', 'tcb', granule, 60);
  var tcgCoef = calcRegression(getSampleImg(img, ref, 'tcg'), ref, 'tcg', 'tcg', granule, 60);
  var tcaCoef = calcRegression(getSampleImg(img, ref, 'tca'), ref, 'tca', 'tca', granule, 60);
  
  return ee.Image(ee.Image.cat(
    applyCoef(img, 'green', greenCoef).toFloat(),
    applyCoef(img, 'red', redCoef).toFloat(),
    applyCoef(img, 'nir', nirCoef).toFloat(),
    applyCoef(img, 'ndvi', nirCoef).toFloat(),
    applyCoef(img, 'tcb', tcbCoef).toFloat(),
    applyCoef(img, 'tcg', tcgCoef).toFloat(),
    applyCoef(img, 'tca', tcaCoef).toFloat())
    .rename(['green', 'red', 'nir', 'ndvi', 'tcb', 'tcg', 'tca'])
    .copyProperties(img, img.propertyNames()));
}



function prepMss(img) {
  var toa = msslib.calcToa(img);
  var toaAddBands = msslib.addTc(msslib.addNdvi(toa));
  var toaAddBandsMask = msslib.applyQaMask(toaAddBands);
  return msslib.applyMsscvm(toaAddBandsMask);
}
exports.prepMss = prepMss;

function processMssWrs1Img(img) {
  var toaAddBandsMsscvmMask = prepMss(img);
  var corrected = correctMssImg(toaAddBandsMsscvmMask);
  return corrected;
}
exports.processMssWrs1Img = processMssWrs1Img;

function processMssWrs1Imgs(params) {
  print('Preparing MSS WRS-1 image processing tasks, please wait');
  var granuleGeom = msslib.getWrs1GranuleGeom(params.wrs1);
  params.aoi = ee.Geometry(granuleGeom.get('centroid'));
  params.wrs = '1';
  var mssCol = msslib.getCol(params)
    .filter(ee.Filter.eq('pr', params.wrs1))
    .map(function(img) {
      return img.set('ref_img', ee.Image(params.baseDir + '/ref'));
    });
  //print(mssCol);
  //var years = ee.List(mssCol.aggregate_array('year').distinct()).sort().getInfo();  // TODO: make this so that all year are written out - include dummies - the script that reads them in assumes that all years exist

  var dummy = ee.Image([0, 0, 0, 0, 0, 0, 0, 0]).selfMask().toShort()
    .rename(['blue', 'green', 'red', 'nir', 'ndvi', 'tcb', 'tcg', 'tca']);
  var outImg;
  for(var y = 1972; y <= 1982; y++) {
    var yrCol = mssCol.filter(ee.Filter.eq('year', y));
    if(yrCol.size().getInfo() === 0) {  // Try to use ee.Algorithms.If - so that the browser does not hang.
      outImg = dummy.set({
          dummy: true,
          year: y,
          'system:time_start': ee.Date.fromYMD(y, 1, 1)
        });
    } else {
      yrCol = yrCol.map(processMssWrs1Img);
      outImg = getMedoid(yrCol, ['green', 'red', 'nir', 'ndvi', 'tcb', 'tcg', 'tca'])
        .set({
          dummy: false,
          year: y,
          'system:time_start': ee.Date.fromYMD(y, 1, 1)
        });
    }
    
    Export.image.toAsset({
      image: outImg,
      description: y.toString(),
      assetId: params.baseDir + '/WRS1_to_WRS2/' + y.toString(),
      region: ee.Feature(granuleGeom.get('granule')).geometry(),
      scale: 60,
      crs: params.crs
    });
  }
}
exports.processMssWrs1Imgs = processMssWrs1Imgs;


function correctMssWrs2(params) {  // NOTE: this is just grabbing 1983 for now.
  var aoi = ee.Feature(
    msslib.getWrs1GranuleGeom(params.wrs1).get('granule')).geometry();
  var mssCol = msslib.getCol({
    aoi: aoi,
    wrs: '2',
    doyRange: params.doyRange
  }).filterDate('1983-01-01', '1984-01-01')
  .map(prepMss);

  return correctMssImgToMedianTm(mssCol, params);
}
exports.correctMssWrs2 = correctMssWrs2;

// #############################################################################
// ### Process TM images ###
// #############################################################################


// Function to get and rename bands of interest from OLI.
function renameOli(img) {
  return img.select(
      ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'pixel_qa'],
      ['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'pixel_qa']);
}

// Function to get and rename bands of interest from ETM+.
function renameTm(img) {
  return img.select(
      ['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'pixel_qa'],
      ['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'pixel_qa']);
}

function tmAddIndices(img) {
  var b = ee.Image(img).select(['blue', 'green', 'red', 'nir', 'swir1', 'swir2']);
  var brt_coeffs = ee.Image.constant([0.2043, 0.4158, 0.5524, 0.5741, 0.3124, 0.2303]);
  var grn_coeffs = ee.Image.constant([-0.1603, -0.2819, -0.4934, 0.7940, -0.0002, -0.1446]);
  var brightness = b.multiply(brt_coeffs).reduce(ee.Reducer.sum()).round().toShort();
  var greenness = b.multiply(grn_coeffs).reduce(ee.Reducer.sum()).round().toShort();
  var angle = (greenness.divide(brightness)).atan().multiply(180 / Math.PI).multiply(100).round().toShort();
  var ndvi = img.normalizedDifference(['nir', 'red']).rename('ndvi').multiply(1000).round().toShort();
  var tc = ee.Image.cat(ndvi, brightness, greenness, angle).rename(['ndvi', 'tcb', 'tcg', 'tca']);
  return img.addBands(tc);
}

function gatherTmCol(params) {
  var granuleGeom = msslib.getWrs1GranuleGeom(params.wrs1);
  var aoi = ee.Feature(granuleGeom.get('granule')).geometry();
  var dateFilter = ee.Filter.calendarRange(params.doyRange[0], params.doyRange[1], 'day_of_year');
  var startDate = ee.Date.fromYMD(params.yearRange[0], 1, 1);
  var endDate = startDate.advance(1, 'year');
  var oliCol = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')
    .filterBounds(aoi).filterDate(startDate, endDate).filter(dateFilter).map(prepOli);
  var etmCol = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')
    .filterBounds(aoi).filterDate(startDate, endDate).filter(dateFilter).map(prepTm);
  var tm5Col = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR')
    .filterBounds(aoi).filterDate(startDate, endDate).filter(dateFilter).map(prepTm);
  var tm4Col = ee.ImageCollection('LANDSAT/LT04/C01/T1_SR')
    .filterBounds(aoi).filterDate(startDate, endDate).filter(dateFilter).map(prepTm);
  return tm4Col.merge(tm5Col).merge(etmCol).merge(oliCol);
}
exports.gatherTmCol = gatherTmCol;

// Define function to prepare OLI images.
function prepOli(img) {
  var orig = img;
  img = renameOli(img);
  img = tmAddIndices(img);
  img = applyCfmask(img);
  return ee.Image(img.copyProperties(orig, orig.propertyNames()));
}

// Define function to prepare ETM+ images.
function prepTm(img) {
  var orig = img;
  img = renameTm(img);
  img = tmAddIndices(img);
  img = applyCfmask(img);
  return ee.Image(img.copyProperties(orig, orig.propertyNames()));
}
exports.prepTm = prepTm;



function getCoincidentTmMssCol(params) {
  var aoi = ee.Feature(
    msslib.getWrs1GranuleGeom(params.wrs1).get('granule')).geometry();
  var mssCol = msslib.getCol({
    aoi: aoi,
    wrs: '2',
    doyRange: params.doyRange
  })
  .map(addTmToMssJoinId);
  
  var tmCol = getTmWrs2Col(aoi).map(addTmToMssJoinId);
  var coincident = coincidentTmMssCol(mssCol, tmCol);
  return coincident;
}
exports.getCoincidentTmMssCol = getCoincidentTmMssCol;

// Function to make normalization function.
function getMss2TmCoefCol(img) {
  var sampleMask = ee.Image.random().gt(0.90);
  var xImg = msslib.addTc(msslib.addNdvi(msslib.calcToa(img)));
  var xImgSamp = xImg.updateMask(sampleMask);
  var yImg = prepTm(ee.Image(xImg.get('coincidentTmMss')));
  var granule= ee.Feature(ee.FeatureCollection('users/jstnbraaten/wrs/wrs2_descending_land')
    .filter(ee.Filter.eq('PR', xImg.getString('pr'))).first()).geometry();
  
  var blueCoef = calcRegression(xImgSamp, yImg, 'green', 'blue', granule, 150); // TODO: this 300 is maybe not ideal - could try a small image sample like 10% or 5% or 1%.
  var greenCoef = calcRegression(xImgSamp, yImg, 'green', 'green', granule, 150);
  var redCoef = calcRegression(xImgSamp, yImg, 'red', 'red', granule, 150);
  var nirCoef = calcRegression(xImgSamp, yImg, 'nir', 'nir', granule, 150);
  var ndviCoef = calcRegression(xImg, yImg, 'ndvi', 'ndvi', granule, 150);
  var tcbCoef = calcRegression(xImg, yImg, 'tcb', 'tcb', granule, 150);
  var tcgCoef = calcRegression(xImg, yImg, 'tcg', 'tcg', granule, 150);
  var tcaCoef = calcRegression(xImg, yImg, 'tca', 'tca', granule, 150);
  
  var coef = {
    blue_coef: blueCoef,
    green_coef: greenCoef,
    red_coef: redCoef,
    nir_coef: nirCoef,
    ndvi_coef: ndviCoef,
    tcb_coef: tcbCoef,
    tcg_coef: tcgCoef,
    tca_coef: tcaCoef,
  };

  xImg = xImg.set('mss_2_tm_coef', coef);
  var xImgCor = _correctMssImg(xImg);
  
  return yImg.select(xImgCor.bandNames()).subtract(xImgCor.select(xImgCor.bandNames()))
    .copyProperties(xImgCor, xImgCor.propertyNames());
}
exports.getMss2TmCoefCol = getMss2TmCoefCol;


function exportTm2mssCoefCol(params) {
  var col = getCoincidentTmMssCol(params);
  var coefFc = col.map(getTm2mssCoefCol);
  Export.table.toAsset({
    collection: coefFc,
    description: 'tm2MssCoefCol',
    assetId: params.baseDir + '/tm2MssCoefCol'
  });
}
exports.exportTm2mssCoefCol = exportTm2mssCoefCol;

function exportMss2TmCoefCol(params) {
  var col = getCoincidentTmMssCol(params);
  var mssOffsetCol = col.map(getMss2TmCoefCol);
  var coefFc = mssOffsetCol.map(function(img) {
    var coefs = ee.Dictionary(img.get('mss_2_tm_coef'));
    var blueCoef = ee.Dictionary(coefs.get('blue_coef'));
    var greenCoef = ee.Dictionary(coefs.get('green_coef'));
    var redCoef = ee.Dictionary(coefs.get('red_coef'));
    var nirCoef = ee.Dictionary(coefs.get('nir_coef'));
    var ndviCoef = ee.Dictionary(coefs.get('ndvi_coef'));
    var tcbCoef = ee.Dictionary(coefs.get('tcb_coef'));
    var tcgCoef = ee.Dictionary(coefs.get('tcg_coef'));
    var tcaCoef = ee.Dictionary(coefs.get('tca_coef'));

    var coef = {
      'blue_slope': blueCoef.getNumber('slope'),
      'blue_intercept': blueCoef.getNumber('intercept'),
      'blue_rmse': blueCoef.getNumber('rmse'),
      'green_slope': greenCoef.getNumber('slope'),
      'green_intercept': greenCoef.getNumber('intercept'),
      'green_rmse': greenCoef.getNumber('rmse'),
      'red_slope': redCoef.getNumber('slope'),
      'red_intercept': redCoef.getNumber('intercept'),
      'red_rmse': redCoef.getNumber('rmse'),
      'nir_slope': nirCoef.getNumber('slope'),
      'nir_intercept': nirCoef.getNumber('intercept'),
      'nir_rmse': nirCoef.getNumber('rmse'),
      'ndvi_slope': ndviCoef.getNumber('slope'),
      'ndvi_intercept': ndviCoef.getNumber('intercept'),
      'ndvi_rmse': ndviCoef.getNumber('rmse'),
      'tcb_slope': tcbCoef.getNumber('slope'),
      'tcb_intercept': tcbCoef.getNumber('intercept'),
      'tcb_rmse': tcbCoef.getNumber('rmse'),
      'tcg_slope': tcgCoef.getNumber('slope'),
      'tcg_intercept': tcgCoef.getNumber('intercept'),
      'tcg_rmse': tcgCoef.getNumber('rmse'),
      'tca_slope': tcaCoef.getNumber('slope'),
      'tca_intercept': tcaCoef.getNumber('intercept'),
      'tca_rmse': tcaCoef.getNumber('rmse'),
    };

    return ee.Feature(ee.Geometry.Point(0, 0)).set(coef)
      .copyProperties(img, ['imgID', 'year', 'path', 'row', 'pr']);
  });
  
  var medianOffset = mssOffsetCol.median().round().toShort();
  var granuleGeom = msslib.getWrs1GranuleGeom(params.wrs1);
  Export.image.toAsset({
    image: medianOffset,
    description: 'medianOffset',
    assetId: params.baseDir + '/mss_offset',
    region: ee.Feature(granuleGeom.get('granule')).geometry(),
    scale: 60,
    crs: params.crs
  });

  Export.table.toAsset({
    collection: coefFc,
    description: 'mss2TmCoefCol',
    assetId: params.baseDir + '/mss2TmCoefCol'
  });
}
exports.exportMss2TmCoefCol = exportMss2TmCoefCol;


function _getMedianCoef(table, coef) {
  return ee.List(table.aggregate_array(coef))
    .reduce(ee.Reducer.median());
}

function getMedianCoef(table) {
  return ee.Dictionary({
    blue_coef: {
      slope: _getMedianCoef(table, 'blue_slope'),
      intercept: _getMedianCoef(table, 'blue_intercept')
    },
    green_coef: {
      slope: _getMedianCoef(table, 'green_slope'),
      intercept: _getMedianCoef(table, 'green_intercept')
    },
    red_coef: {
      slope: _getMedianCoef(table, 'red_slope'),
      intercept: _getMedianCoef(table, 'red_intercept')
    },
    nir_coef: {
      slope: _getMedianCoef(table, 'nir_slope'),
      intercept: _getMedianCoef(table, 'nir_intercept')
    },
    ndvi_coef: {
      slope: _getMedianCoef(table, 'ndvi_slope'),
      intercept: _getMedianCoef(table, 'ndvi_intercept')
    },
    tcb_coef: {
      slope: _getMedianCoef(table, 'tcb_slope'),
      intercept: _getMedianCoef(table, 'tcb_intercept')
    },
    tcg_coef: {
      slope: _getMedianCoef(table, 'tcg_slope'),
      intercept: _getMedianCoef(table, 'tcg_intercept')
    },
    tca_coef: {
      slope: _getMedianCoef(table, 'tca_slope'),
      intercept: _getMedianCoef(table, 'tca_intercept')
    }
  });
}
exports.getMedianCoef = getMedianCoef;

function _correctMssImg(img) {
  var coefs = ee.Dictionary(img.get('mss_2_tm_coef'));
  return ee.Image(ee.Image.cat(
    applyCoef(img, 'green', ee.Dictionary(coefs.get('blue_coef'))).toFloat(),
    applyCoef(img, 'green', ee.Dictionary(coefs.get('green_coef'))).toFloat(),
    applyCoef(img, 'red', ee.Dictionary(coefs.get('red_coef'))).toFloat(),
    applyCoef(img, 'nir', ee.Dictionary(coefs.get('nir_coef'))).toFloat(),
    applyCoef(img, 'ndvi', ee.Dictionary(coefs.get('ndvi_coef'))).toFloat(),
    applyCoef(img, 'tcb', ee.Dictionary(coefs.get('tcb_coef'))).toFloat(),
    applyCoef(img, 'tcg', ee.Dictionary(coefs.get('tcg_coef'))).toFloat(),
    applyCoef(img, 'tca', ee.Dictionary(coefs.get('tca_coef'))).toFloat())
    .rename(['blue', 'green', 'red', 'nir', 'ndvi', 'tcb', 'tcg', 'tca'])
    .copyProperties(img, img.propertyNames()));
}
exports._correctMssImg = _correctMssImg;

function correctMssImgToMedianTm(col, params) {
  var table = ee.FeatureCollection(params.baseDir + '/mss2TmCoefCol');
  var offset = ee.Image(params.baseDir + '/mss_offset');
  var coefs = getMedianCoef(table);
  return col.map(function(img) {
    return img
      .set('mss_2_tm_coef', coefs);
  })
  .map(_correctMssImg)
  .map(function(img) {
    return img.add(offset).round().toShort().copyProperties(img, img.propertyNames());  // NOTE: no offset - img.round().toShort().copyProperties(img, img.propertyNames());
  });
}
exports.correctMssImgToMedianTm = correctMssImgToMedianTm;

function getFinalCorrectedMssCol(params) {
  var mssCol = ee.ImageCollection([]);
  for(var y = 1972; y <= 1982; y++) {
    var img = ee.Image(params.baseDir + '/WRS1_to_WRS2/' + y.toString())
      .set('system:time_start', ee.Date.fromYMD(y, 1 ,1).millis());
    mssCol = mssCol.merge(ee.ImageCollection(img));
  }
  return correctMssImgToMedianTm(mssCol, params);
}

function exportFinalCorrectedMssCol(params) {
  var mssCol = getFinalCorrectedMssCol(params);
  var granuleGeom = msslib.getWrs1GranuleGeom(params.wrs1);
  // var dummy = ee.Image([0, 0, 0, 0, 0, 0, 0, 0]).selfMask().toShort()
  //   .rename(['blue', 'green', 'red', 'nir', 'ndvi', 'tcb', 'tcg', 'tca']);

  for(var y = 1972; y <= 1982; y++) {
    var yrCol = mssCol.filter(ee.Filter.eq('year', y));
    
    // // Deal with missing years - provide a dummy.
    // var outImg = ee.Algorithms.If({
    //   condition: yrCol.size(),
    //   trueCase: yrCol.first().resample('bicubic'),  // NOTE: not sure about the resample?,
    //   falseCase: dummy.set('system:time_start', ee.Date.fromYMD(y, 1 ,1).millis())
    // });

    Export.image.toAsset({
      image: yrCol.first().resample('bicubic'), //outImg,
      description: y.toString(),
      assetId: params.baseDir + '/WRS1_to_TM/' + y.toString(),
      region: ee.Feature(granuleGeom.get('granule')).geometry(),
      scale: 30,
      crs: params.crs,
      maxPixels: 1e13
    });
  }
}
exports.exportFinalCorrectedMssCol = exportFinalCorrectedMssCol;


// #############################################################################
// ### Final collection assembly ###
// #############################################################################

function getColForLandTrendrOnTheFly(params) { // Does not rely on WRS1_to_TM assets
  var mssCol = getFinalCorrectedMssCol(params);

  var tmCol = ee.ImageCollection([]);
  for(var y=1983; y<=2020; y++) {
    params.yearRange = [y, y];
    var thisYearCol = getMedoid(gatherTmCol(params), ['blue', 'green', 'red', 'nir', 'ndvi', 'tcb', 'tcg', 'tca'])
      .set('system:time_start', ee.Date.fromYMD(y, 1 ,1).millis());
    tmCol = tmCol.merge(ee.ImageCollection(thisYearCol.toShort()));
  }
  
  var combinedCol = mssCol.merge(tmCol).map(function(img) {
    return img.select('ndvi').multiply(-1).rename('LTndvi')  // TODO: move this into the run landtrandr function - get the fitting index from params
      .addBands(img.select(['blue', 'green', 'red', 'nir', 'ndvi', 'tcb', 'tcg', 'tca']))  // TODO: move this into the run landtrandr function - what indices should be FTV
      .set('system:time_start', img.get('system:time_start'));
  }).sort('system:time_start');

  return combinedCol;
}
exports.getColForLandTrendrOnTheFly = getColForLandTrendrOnTheFly;

function getColForLandTrendrFromAsset(params) { // Relies on WRS1_to_TM assets
  var mssCol = ee.ImageCollection([]);
  for(var y=1972; y<=1982; y++) {
    var img = ee.Image(params.baseDir + '/WRS1_to_TM/' + y.toString())
      .set('system:time_start', ee.Date.fromYMD(y, 1 ,1).millis());
    mssCol = mssCol.merge(ee.ImageCollection(img));
  }
  
  var mss1983 = ee.ImageCollection(getMedoid(correctMssWrs2(params), ['blue', 'green', 'red', 'nir', 'ndvi', 'tcb', 'tcg', 'tca'])
    .set('system:time_start', ee.Date.fromYMD(1983, 1 ,1).millis()));
  
  var tmCol = ee.ImageCollection([]);
  for(var y=1984; y<=2020; y++) {
    params.yearRange = [y, y];
    var thisYearCol = getMedoid(gatherTmCol(params), ['blue', 'green', 'red', 'nir', 'ndvi', 'tcb', 'tcg', 'tca'])
      .set('system:time_start', ee.Date.fromYMD(y, 1 ,1).millis());
    tmCol = tmCol.merge(ee.ImageCollection(thisYearCol.toShort()));
  }
  
  var combinedCol = mssCol.merge(mss1983).merge(tmCol).map(function(img) {
    return img.select('ndvi').multiply(-1).rename('LTndvi')  // TODO: move this into the run landtrandr function - get the fitting index from params
      .addBands(img.select(['blue', 'green', 'red', 'nir', 'ndvi', 'tcb', 'tcg', 'tca']))  // TODO: move this into the run landtrandr function - what indices should be FTV
      .set('system:time_start', img.get('system:time_start'));
  }).sort('system:time_start');

  return combinedCol;
}
exports.getColForLandTrendrFromAsset = getColForLandTrendrFromAsset;


function runLandTrendrMss2Tm(params) {
  var ltCol = getColForLandTrendrFromAsset(params); // alternative: getColForLandTrendrOnTheFly(params)
  var lt = ee.Algorithms.TemporalSegmentation.LandTrendr({
    timeSeries: ltCol,
    maxSegments: 10,
    spikeThreshold: 0.7,
    vertexCountOvershoot: 3,
    preventOneYearRecovery: true,
    recoveryThreshold: 0.5,
    pvalThreshold: 0.05,
    bestModelProportion: 0.75,
    minObservationsNeeded: 6
  });
  return lt;
}
exports.runLandTrendrMss2Tm = runLandTrendrMss2Tm;

// #############################################################################
// ### Functions under development (Annie Taylor) ###
// #############################################################################

function displayCollection(col) {
  var rgbviz = {
    bands: ['red','green','blue'],
    min: 100,
    max: 2000,
    gamma: [1.2]
  };
  Map.centerObject(col.first(), 8);
  Map.addLayer(col, rgbviz, 'Full Landsat Collection',false);
}
exports.displayCollection = displayCollection; 


function animateCollection(col) {
  var rgbviz = {
    bands: ['red','green','blue'],
    min: 100,
    max: 2000,
    gamma: [1.2]
  };
  // TODO: add year of image as label in animation
  // col = col.map(function(img) {
  //   img = img.set({label: ee.String(img.get('system:id'))})
  //   return img
  // })
  Map.centerObject(col.first(), 8);
  // run the animation
  animation.animate(col, {
    vis: rgbviz,
    timeStep: 1500,
    maxFrames: col.size()
  })
}
exports.animateCollection = animateCollection; 

function displayGreatestDisturbance(lt, params) {
  var granuleGeom = ee.Feature(msslib.getWrs1GranuleGeom(params.wrs1)
    .get('granule')).geometry();
  
  var currentYear = new Date().getFullYear();  // TODO: make sure there is not a better way to get year from image metadata eg
  var changeParams = { // TODO: allow a person to override these params
    delta:  'loss',
    sort:   'greatest',
    year:   {checked:true, start:1972, end:currentYear},  // TODO: make sure there is not a better way to get years from image metadata eg
    mag:    {checked:true, value:200,  operator:'>'},
    dur:    {checked:true, value:4,    operator:'<'},
    preval: {checked:true, value:300,  operator:'>'},
    mmu:    {checked:true, value:11},
  };
  // Note: add index to changeParams object this is hard coded to NDVI because currently that is the only option.
  changeParams.index = 'NDVI';
  var changeImg = ltgee.getChangeMap(lt, changeParams);
  var palette = ['#9400D3', '#4B0082', '#0000FF', '#00FF00',
                  '#FFFF00', '#FF7F00', '#FF0000'];
  var yodVizParms = {
    min: 1972, // TODO: make sure there is not a better way to get year from image metadata eg
    max: currentYear, // TODO: make sure there is not a better way to get year from image metadata eg
    palette: palette
  };
  var magVizParms = {
    min: 200,
    max: 800,
    palette: palette
  };
  Map.centerObject(granuleGeom, 12);  // Zoom in pretty far otherwise the mmu filter is going to take forever (probably crash)
  // display two change attributes to map
  Map.addLayer(changeImg.select(['mag']), magVizParms, 'Magnitude of Change');
  Map.addLayer(changeImg.select(['yod']), yodVizParms, 'Year of Detection');
}
exports.displayGreatestDisturbance = displayGreatestDisturbance;

// #############################################################################
// ### TM to MSS functions ###
// #############################################################################



// // Function to make normalization function.
// function getTm2mssCoefCol(img) { //function makeCorrectionFun(refImgPath) {
//   //var img = coCol.first();
//   var yImg = msslib.addNdvi(msslib.calcToa(img));
//   var xImg = ee.Image(yImg.get('coincidentTmMss'));
  
//   var xImgNdvi = xImg.normalizedDifference(['B4', 'B3']).rename(['ndvi']);
//   var mask = getCfmask(xImg);
//   xImg = xImg.addBands(xImgNdvi).updateMask(mask);
  
//   var greenCoef = calcRegression(xImg, yImg, 'B2', 'green');
//   var redCoef = calcRegression(xImg, yImg, 'B3', 'red');
//   var nirCoef = calcRegression(xImg, yImg, 'B4', 'nir');
//   var ndviCoef = calcRegression(xImg, yImg, 'ndvi', 'ndvi');
  
//   var granuleGeoms = msslib.getWrs1GranuleGeom('049029');  // TODO - NEED TO GET THIS FROM THE PARAMS - could really just make this a dummy geom at 0, 0 - it's not needed.
//   var centroid = ee.Geometry(granuleGeoms.get('centroid'));
  
//   return ee.Feature(centroid).set(
//     {
//       'green_slope': greenCoef.getNumber('slope'),
//       'green_intercept': greenCoef.getNumber('intercept'),
//       'red_slope': redCoef.getNumber('slope'),
//       'red_intercept': redCoef.getNumber('intercept'),
//       'nir_slope': nirCoef.getNumber('slope'),
//       'nir_intercept': nirCoef.getNumber('intercept'),
//       'ndvi_slope': ndviCoef.getNumber('slope'),
//       'ndvi_intercept': ndviCoef.getNumber('intercept')
//     }).copyProperties(yImg, ['imgID', 'year', 'path', 'row', 'pr']); // TODO: need to add RMSE. - Could just exclude the system:footprint
// }
// exports.getTm2mssCoefCol = getTm2mssCoefCol;


// function _correctTmImg(img) {
//   var coefs = ee.Dictionary(img.get('tm_2_mss_coef'));
//   return ee.Image(ee.Image.cat(
//     applyCoef(img, 'B2', ee.Dictionary(coefs.get('green_coef'))).toFloat(),
//     applyCoef(img, 'B3', ee.Dictionary(coefs.get('red_coef'))).toFloat(),
//     applyCoef(img, 'B4', ee.Dictionary(coefs.get('nir_coef'))).toFloat(),
//     applyCoef(img, 'ndvi', ee.Dictionary(coefs.get('ndvi_coef'))).toFloat())
//     .rename(['green', 'red', 'nir', 'ndvi'])
//     .copyProperties(img, img.propertyNames()));
// }

// function correctTmImg(col, params) {
//   var table = ee.FeatureCollection(params.baseDir + '/tm2MssCoefCol');
//   var coefs = getMedianCoef(table);
//   return col.map(function(img) {
//     return img.set('tm_2_mss_coef', coefs);
//   })
//   .map(_correctTmImg);
// }
// exports.correctTmImg = correctTmImg;


// function runLandTrendrTm2Mss() {
//   var mssCol = ee.ImageCollection([]);
//   for(var y=1972; y<=1982; y++) {
//     var img = ee.Image("users/braaten/llr/test_proj_1/049029/" + y.toString())
//       .set('system:time_start', ee.Date.fromYMD(y, 1 ,1).millis());
//     mssCol = mssCol.merge(ee.ImageCollection(img));
//   }

//   var tmCol = ee.ImageCollection([]);
//   for(var y=1983; y<=2020; y++) {
//     params.yearRange = [y, y];
//     var thisYear = gatherTmCol(params);
//     var thisYearCol = correctTmImg(thisYear, params).median()  // TODO: this needs to be medoid.
//       .set('system:time_start', ee.Date.fromYMD(y, 1 ,1).millis());
//     tmCol = tmCol.merge(ee.ImageCollection(thisYearCol));
//   }
  
//   var ltCol = mssCol.merge(tmCol).map(function(img) {
//     return img.select('ndvi').multiply(-1).rename('LTndvi')
//       .addBands(img.select(['green', 'red', 'nir']))
//       .set('system:time_start', img.get('system:time_start'));
//   });
  
//   var lt = ee.Algorithms.TemporalSegmentation.LandTrendr({
//     timeSeries: ltCol,
//     maxSegments: 10,
//     spikeThreshold: 0.7,
//     vertexCountOvershoot: 3,
//     preventOneYearRecovery: true,
//     recoveryThreshold: 0.5,
//     pvalThreshold: 0.05,
//     bestModelProportion: 0.75,
//     minObservationsNeeded: 6
//   });

//   return lt;
// }
