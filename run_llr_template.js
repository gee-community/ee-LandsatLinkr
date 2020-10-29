// #############################################################################
// #############################################################################
// #############################################################################

/**
 * 1. View WRS-1 granules - figure out what WRS-1 granule to process
 * -- Make a processing dir: https://gist.github.com/jdbcode/36f5a04329d5d85c43c0408176c51e6d
 * 2. Create MSS WRS-1 reference image - for MSS WRS1 to MSS WRS2 harmonization
 * 3. View WRS-1 collection - identify bad MSS images
 * 4. Prepare MSS WRS-1 images
 * 5. Get TM-to-MSS correction coefficients
 * 6. Export MSS-to-TM corrected images
 * 7. Inspect the full time series collection - builds col - add code lines to add col to Map and use inspector etc to explore time series
 * 8. Run LandTrendr - returns the raw output from LandTrendr - you'll need to add lines to run lt-gee functions or export e.g.
 */

var LLR_STEP = 1;

// #############################################################################

var PROJ_PATH = 'users/braaten/LandTrendr';   // Must be the same path used to create the asset folder - cannot contain / at end - check for this in the code.
var WRS_1_GRANULE = '049030';
var CRS = 'EPSG:3857';

var DOY_RANGE = [160, 254];
var MAX_CLOUD = 50;
var MAX_GEOM_RMSE = 0.5;

var EXCLUDE_IDS = [
  'LM10490301972210AAA05',
  'LM20490301975167AAA02',
  'LM20490301975185AAA02',
  'LM20490301976216GDS01',
  'LM20490301976252GDS01',
  'LM30490301978196GDS03',
  'LM20490301978241AAA02',
  'LM20490301979200AAA05',
  'LM20490301981171AAA03'
];

// #############################################################################
// #############################################################################
// #############################################################################

var params = {
  maxRmseVerify: MAX_GEOM_RMSE,
  maxCloudCover: MAX_CLOUD,
  doyRange: DOY_RANGE,
  wrs1: WRS_1_GRANULE,
  crs: CRS,
  excludeIds: EXCLUDE_IDS,
  baseDir: PROJ_PATH + '/' + WRS_1_GRANULE
};

var llr = require('users/jstnbraaten/modules:landsatlinkr/landsatlinkr.js');
switch (LLR_STEP) {
  case 1:
    llr.wrs1GranuleSelector();
    break;
  case 2:
    llr.exportMssRefImg(params);
    break;
  case 3:
    llr.viewWrs1Col(params);
    break;
  case 4:
    llr.processMssWrs1Imgs(params);
    break;
  case 5:
    llr.exportMss2TmCoefCol(params);
    break;
  case 6:
    llr.exportFinalCorrectedMssCol(params);
    break;
  case 7:
    var col = llr.getColForLandTrendrFromAsset(params);
    break;
  case 8:
    var lt = llr.runLandTrendrMss2Tm(params);
    break;
}
