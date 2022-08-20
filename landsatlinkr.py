from datetime import date
import copy
import os
import time
import math
import pprint
from IPython.display import Image
import subprocess
import sys
import ee

from traitlets.traitlets import default
def getPr(img):
    path = img.getNumber('WRS_PATH').format('%03d')
    row = img.getNumber('WRS_ROW').format('%03d')
    return path.cat(row)


def filterById_doit(id, col):
  return ee.ImageCollection(col).filter(
      ee.Filter.neq('LANDSAT_SCENE_ID', ee.String(id)))


def filterById(col, imgList):
  return ee.ImageCollection(ee.List(imgList).iterate(filterById_doit, col))


def filterCol(col, params, wrs):
  # Adjust band present property names depending on WRS (1 or 2).
  bandsPresent = {
    'wrs1': [
      'PRESENT_BAND_4', 'PRESENT_BAND_5', 'PRESENT_BAND_6', 'PRESENT_BAND_7'
    ],
    'wrs2': [
      'PRESENT_BAND_1', 'PRESENT_BAND_2', 'PRESENT_BAND_3', 'PRESENT_BAND_4'
    ]
  }

  if params['aoi']:
    col = col.filterBounds(params['aoi'])

  col = (col.filter(ee.Filter.neq('DATA_TYPE', 'L1G'))
            .filter(ee.Filter.eq(bandsPresent[wrs][0], 'Y'))
            .filter(ee.Filter.eq(bandsPresent[wrs][1], 'Y'))
            .filter(ee.Filter.eq(bandsPresent[wrs][2], 'Y'))
            .filter(ee.Filter.eq(bandsPresent[wrs][3], 'Y'))
            .filter(ee.Filter.lte('GEOMETRIC_RMSE_VERIFY', params['maxRmseVerify']))
            .filter(ee.Filter.lte('CLOUD_COVER', params['maxCloudCover'])))

  if params['yearRange']:
    col = col.filter(ee.Filter.calendarRange(
        params['yearRange'][0], params['yearRange'][1], 'year'))

  if params['doyRange']:
    col = col.filter(ee.Filter.calendarRange(
        params['doyRange'][0], params['doyRange'][1], 'day_of_year'))

  if params['excludeIds']:
    col = filterById(col, params['excludeIds'])

  return col


def getCol(params):
  # Define default filter parameters.
  _params = {
    'aoi': None,
    'maxRmseVerify': 0.5,
    'maxCloudCover': 50,
    'wrs': '1&2',
    'yearRange': [1972, 2000],
    'doyRange': [1, 365],
    'excludeIds': None
  }

  # Replace default params with provided params.
  if params:
    for param in params:
      _params[param] = params[param] or _params[param]

  # Initialize WRS-1 and WRS-2 collections.
  wrs1Col = ee.ImageCollection([])
  wrs2Col = ee.ImageCollection([])

  # Gather MSS WRS-1 images, filter as requested, designate as 'WRS-1'.
  if '1' in _params['wrs']:
    mss1T1 = filterCol(
        ee.ImageCollection('LANDSAT/LM01/C02/T1'), _params, 'wrs1')
    mss1T2 = filterCol(
        ee.ImageCollection('LANDSAT/LM01/C02/T2'), _params, 'wrs1')
    mss2T1 = filterCol(
        ee.ImageCollection('LANDSAT/LM02/C02/T1'), _params, 'wrs1')
    mss2T2 = filterCol(
        ee.ImageCollection('LANDSAT/LM02/C02/T2'), _params, 'wrs1')
    mss3T1 = filterCol(
        ee.ImageCollection('LANDSAT/LM03/C02/T1'), _params, 'wrs1')
    mss3T2 = filterCol(
        ee.ImageCollection('LANDSAT/LM03/C02/T2'), _params, 'wrs1')
    wrs1Col = (ee.ImageCollection(ee.FeatureCollection(
        [mss1T1, mss1T2, mss2T1, mss2T2, mss3T1, mss3T2]).flatten())
        .select(['B.|QA_PIXEL|QA_RADSAT'],
                ['green', 'red', 'red-edge', 'nir', 'QA_PIXEL', 'QA_RADSAT'])
        .map(lambda img: img.set('wrs', 'WRS-1')))

  # Gather MSS WRS-2 images, filter as requested, designate as 'WRS-2'.
  if '2' in _params['wrs']:
    mss4T1 = filterCol(
        ee.ImageCollection('LANDSAT/LM04/C02/T1'), _params, 'wrs2');
    mss4T2 = filterCol(
        ee.ImageCollection('LANDSAT/LM04/C02/T2'), _params, 'wrs2');
    mss5T1 = filterCol(
        ee.ImageCollection('LANDSAT/LM05/C02/T1'), _params, 'wrs2');
    mss5T2 = filterCol(
        ee.ImageCollection('LANDSAT/LM05/C02/T2'), _params, 'wrs2');
    wrs2Col = (ee.ImageCollection(ee.FeatureCollection(
        [mss4T1, mss4T2, mss5T1, mss5T2]).flatten())
        .select(['B.|QA_PIXEL|QA_RADSAT'],
                ['green', 'red', 'red-edge', 'nir', 'QA_PIXEL', 'QA_RADSAT'])
        .map(lambda img: img.set('wrs', 'WRS-2')))

  # Return time-sorted, merged, WRS-1 and WRS-2 collection with filter params
  # attached.
  return ee.ImageCollection(ee.FeatureCollection([wrs1Col, wrs2Col]).flatten()).map(lambda img: img.set({
      'start_doy': _params['doyRange'][0],
      'end_doy': _params['doyRange'][1],
      'year': img.date().get('year'),
      'doy': img.date().getRelative('day', 'year'),
      'pr': getPr(img)
      # composite_year:  # TODO
    })).sort('system:time_start')


def getWrs1GranuleGeom(granuleId):
    granule = ee.Feature(
            ee.FeatureCollection('users/jstnbraaten/wrs/wrs1_descending_land')
            .filter(ee.Filter.eq('PR', granuleId)).first())
    centroid = granule.centroid(300).geometry(300)
    bounds = granule.geometry(300).buffer(40000)
    return ee.Dictionary({
        'granule': granule,
        'centroid': centroid,
        'bounds': bounds
    })


def scaleDn(img, unit):
  mult = 'REFLECTANCE_MULT_BAND'
  add = 'REFLECTANCE_ADD_BAND'
  if unit == 'radiance':
    mult = 'RADIANCE_MULT_BAND'
    add = 'RADIANCE_ADD_BAND'

  gainBands = (ee.List(img.propertyNames())
                      .filter(ee.Filter.stringContains('item', mult))
                      .sort())
  biasBands = (ee.List(img.propertyNames())
                      .filter(ee.Filter.stringContains('item', add))
                      .sort())

  gainImg = ee.Image.cat(
      ee.Image.constant(img.get(gainBands.getString(0))),
      ee.Image.constant(img.get(gainBands.getString(1))),
      ee.Image.constant(img.get(gainBands.getString(2))),
      ee.Image.constant(img.get(gainBands.getString(3)))).toFloat()

  biasImg = ee.Image.cat(
      ee.Image.constant(img.get(biasBands.getString(0))),
      ee.Image.constant(img.get(biasBands.getString(1))),
      ee.Image.constant(img.get(biasBands.getString(2))),
      ee.Image.constant(img.get(biasBands.getString(3)))).toFloat()

  dnImg = img.select([0, 1, 2, 3]).multiply(gainImg).add(biasImg).toFloat()

  return img.addBands(dnImg, None, True)


def calcToa(img):
  return scaleDn(img, 'reflectance')


def addTc(img):
  bands = img.select([0, 1, 2, 3])
  tcbCoeffs = ee.Image.constant([0.433, 0.632, 0.586, 0.264])
  tcgCoeffs = ee.Image.constant([-0.290, -0.562, 0.600, 0.491])
  tcyCoeffs = ee.Image.constant([-0.829, 0.522, -0.039, 0.194])
  tcb = bands.multiply(tcbCoeffs).reduce(ee.Reducer.sum()).toFloat()
  tcg = bands.multiply(tcgCoeffs).reduce(ee.Reducer.sum()).toFloat()
  tcy = bands.multiply(tcyCoeffs).reduce(ee.Reducer.sum()).toFloat()
  tca = (tcg.divide(tcb)).atan().multiply(180 / math.pi).toFloat()
  tc = ee.Image.cat(tcb, tcg, tcy, tca).rename('tcb', 'tcg', 'tcy', 'tca')
  return ee.Image(img.addBands(tc).copyProperties(img, img.propertyNames()))


def addNdvi(img):
  ndvi = img.normalizedDifference(['nir', 'red']).rename('ndvi')
  return ee.Image(img.addBands(ndvi).copyProperties(img, img.propertyNames()))


def scaleDn(img, unit):
  mult = 'REFLECTANCE_MULT_BAND'
  add = 'REFLECTANCE_ADD_BAND'
  if unit == 'radiance':
    mult = 'RADIANCE_MULT_BAND'
    add = 'RADIANCE_ADD_BAND'

  gainBands = (ee.List(img.propertyNames())
                      .filter(ee.Filter.stringContains('item', mult))
                      .sort())
  biasBands = (ee.List(img.propertyNames())
                      .filter(ee.Filter.stringContains('item', add))
                      .sort())

  gainImg = ee.Image.cat(
      ee.Image.constant(img.get(gainBands.getString(0))),
      ee.Image.constant(img.get(gainBands.getString(1))),
      ee.Image.constant(img.get(gainBands.getString(2))),
      ee.Image.constant(img.get(gainBands.getString(3)))).toFloat()

  biasImg = ee.Image.cat(
      ee.Image.constant(img.get(biasBands.getString(0))),
      ee.Image.constant(img.get(biasBands.getString(1))),
      ee.Image.constant(img.get(biasBands.getString(2))),
      ee.Image.constant(img.get(biasBands.getString(3)))).toFloat()

  dnImg = img.select([0, 1, 2, 3]).multiply(gainImg).add(biasImg).toFloat()

  return img.addBands(dnImg, None, True)


def calcRad(img):
  return scaleDn(img, 'radiance')


def calcToa(img):
  return scaleDn(img, 'reflectance')


visDn = {
  'bands': ['nir', 'red', 'green'], 
  'min': [47, 20, 27],
  'max': [142, 92, 71],
  'gamma': [1.2, 1.2, 1.2]
}


visRad = {
  'bands': ['nir', 'red', 'green'],
  'min': [23, 15, 25],
  'max': [67, 62, 64],
  'gamma': [1.2, 1.2, 1.2]
}


visToa = {
  'bands': ['nir', 'red', 'green'],
  'min': [0.0896, 0.0322, 0.0464],
  'max': [0.2627, 0.1335, 0.1177],
  'gamma': [1.2, 1.2, 1.2]
}


visNdvi = {
  'bands': ['ndvi'], 'min': 0.1, 'max': 0.8
}

def viewThumbnails(col, params=None):
  print('Please wait patiently, images may not load immediately\n')

  _params = {
    'unit': 'toa',
    'display': 'nir|red|green',
    'visParams': None
  }

  if params:
    for param in params:
      _params[param] = params[param] or _params[param]
  
  settings = {
    'unit': {
      'dn': lambda img: img,
      'rad': calcRad,
      'toa': calcToa
    },
    'display': {
      'nir|red|green': {
        'dn': visDn,
        'rad': visRad,
        'toa': visToa  
      },
      'ndvi': {
        'dn': visNdvi,
        'rad': visNdvi,
        'toa': visNdvi
      }
    }
  }

  nImgs = col.size().getInfo()
  imgList = col.sort('system:time_start').toList(nImgs).getInfo()
  for i in range(0, nImgs):
    id = imgList[i]['id']
    img = applyQaMask(ee.Image(id)).select(['B.'], ['green', 'red', 'red-edge', 'nir'])
    img = settings['unit'][_params['unit']](img)
    if _params['display'] == 'ndvi':
      img = addNdvi(img)

    visParams = settings['display'][_params['display']][_params['unit']]
    if _params['visParams']:
      visParams = _params['visParams']
    
    imgVis = img.unmask(0).visualize(**visParams)
    date = img.date().format('YYYY-MM-dd').getInfo()
    sceneId = img.get('LANDSAT_SCENE_ID').getInfo()
    print(f'Date: {date} | Scene ID: {sceneId}')
    display(Image(url=imgVis.getThumbURL({
        'dimensions': 512,
        'crs': 'EPSG:3857'})))
    print('\n')


def getQaMask(img):
  qaPixelMask = img.select('QA_PIXEL').bitwiseAnd(int('11111', 2)).eq(0)
  qaRadsatMask = img.select('QA_RADSAT').eq(0)
  return qaPixelMask.updateMask(qaRadsatMask).rename('QA_mask')


def applyQaMask(img): # TODO: I don't think this is being applied except for thumbVis
  return img.updateMask(getQaMask(img))


def getDem(img):
    aw3d30 = ee.Image('JAXA/ALOS/AW3D30/V2_2').select('AVE_DSM').rename('elev')
    GMTED2010 = ee.Image('USGS/GMTED2010').rename('elev')
    return ee.ImageCollection([GMTED2010, aw3d30]) \
      .mosaic() \
      .reproject(img.projection())

def waterLayer(img):
    # Threshold on NDVI.
    mssWater = img.normalizedDifference(['nir', 'red']).lt(-0.085)

    # Get max extent of water 1985-2018.
    waterExtent = ee.Image('JRC/GSW1_1/GlobalSurfaceWater').select('max_extent')

    # Get intersection of MSS water and max extent.
    return mssWater.multiply(waterExtent) \
      .reproject(img.projection()) \
      .rename('water')

def radians(img):
    return img.toFloat().multiply(math.pi).divide(180)

def getIll(img, slope, aspect):
    # Get sun info.
    azimuth = img.get('SUN_AZIMUTH')
    zenith = ee.Number(90).subtract(img.getNumber('SUN_ELEVATION'))

    # Convert slope and aspect degrees to radians.
    slopeRad = radians(slope)
    aspectRad = radians(aspect)

    # Calculate illumination.
    azimuthImg = radians(ee.Image.constant(azimuth))
    zenithImg = radians(ee.Image.constant(zenith))
    left = zenithImg.cos().multiply(slopeRad.cos())
    right = zenithImg.sin() \
                  .multiply(slopeRad.sin()) \
                  .multiply(azimuthImg.subtract(aspectRad).cos())
    return left.add(right)

def topoCorrB4(img, dem):
    # Get terrain layers.
    terrain = ee.Algorithms.Terrain(dem)
    slope = terrain.select(['slope'])
    aspect = terrain.select(['aspect'])

    # Get k image.
    # define polynomial coefficients to calc Minnaert value as function of slope
    # Ge, H., Lu, D., He, S., Xu, A., Zhou, G., & Du, H. (2008). Pixel-based
    # Minnaert correction method for reducing topographic effects on a Landsat 7
    # ETM+ image. Photogrammetric Engineering & Remote Sensing, 74(11),
    # 1343-1350. |
    # https:#orst.library.ingentaconnect.com/content/asprs/pers/2008/00000074/00000011/art00003?crawler=True&mimetype=application/pdf
    kImg = (slope.resample('bilinear')
                 .where(
                     slope.gt(50),
                     50)  # Set max slope at 50 degrees - paper does not sample \
                 .polynomial([
                   1.0021313684, -0.1308793751, 0.0106861276, -0.0004051135,
                   0.0000071825, -4.88e-8
                 ]))

    # Get illumination.
    ill = getIll(img, slope, aspect)

    # Correct NIR reflectance for topography.
    cosTheta = radians(ee.Image.constant(ee.Number(90).subtract(
                             ee.Number(img.get('SUN_ELEVATION'))))).cos()
    correction = (cosTheta.divide(ill)).pow(kImg)
    return img.select('nir').multiply(correction)

def cloudLayer(img):
    # Identify cloud pixels.
    cloudPixels = (img.normalizedDifference(['green', 'red'])
                        .gt(0)
                        .multiply(img.select('green').gt(0.175))
                        .add(img.select('green').gt(0.39))
                        .gt(0))

    # Nine-pixel minimum connected component sieve.
    cloudPixels = (cloudPixels.selfMask()
                    .connectedPixelCount(10, True)
                    .reproject(img.projection())
                    .gte(0)
                    .unmask(0)
                    .rename('cloudtest'))

    # Define kernel for buffer.
    kernel = ee.Kernel.circle(**{'radius': 2, 'units': 'pixels', 'normalize': True})

    # Two pixel buffer, eight neighbor rule.
    return (cloudPixels.focalMax(**{'radius': 2, 'kernel': kernel})
      .reproject(img.projection())
      .rename('clouds'))

def shadowLayer(img, dem, clouds):
    # Correct B4 reflectance for topography.
    b4c = topoCorrB4(img, dem)

    # Threshold B4 - target dark pixels.
    shadows = b4c.lt(0.11);  # Make this True for all pixels to use full cloud projection.

    # Project clouds as potential shadow.
    shadow_azimuth = ee.Number(90).subtract(ee.Number(img.get('SUN_AZIMUTH')))
    cloudProj = (clouds.directionalDistanceTransform(shadow_azimuth, 50)
                      .reproject(**{'crs': img.projection(), 'scale': 60})
                      .select('distance')
                      .gt(0)
                      .unmask(0))

    # Get water layer.
    water = waterLayer(img)

    # Exclude water pixels from intersection of cloud projection and dark pixels.
    return (shadows.multiply(water.Not())
      .multiply(cloudProj)
      .focalMax(2)
      .reproject(img.projection()))

def applyMsscvm(img):
    dem = getDem(img)
    water = waterLayer(img)
    b4c = topoCorrB4(img, dem)
    clouds = cloudLayer(img)
    shadows = shadowLayer(img, dem, clouds)
    mask = clouds.add(shadows).eq(0)
    return img.updateMask(mask)


msslib = {
    'getWrs1GranuleGeom': getWrs1GranuleGeom,
    'getCol': getCol,
    'calcToa': calcToa,
    'addTc': addTc,
    'addNdvi': addNdvi,
    'viewThumbnails': viewThumbnails,
    'applyQaMask': applyQaMask,
    'applyMsscvm': applyMsscvm
}

params = None #  GLOBAL dict redefined later.

def getTmWrs2Col(aoi):
    tm4 = ee.ImageCollection("LANDSAT/LT04/C02/T1_L2") \
        .filterBounds(aoi)
    tm5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2") \
        .filterBounds(aoi)
    return tm4.merge(tm5)

def coincidentTmMssCol(mssWrs2Col, tmWrs2Col):
    filter = ee.Filter.equals(**{'leftField': 'imgID', 'rightField': 'imgID'})
    join = ee.Join.saveFirst('coincidentTmMss')
    return ee.ImageCollection(join.apply(mssWrs2Col, tmWrs2Col, filter))

def addTmToMssJoinId(img):
    date = ee.Image(img).date()
    year = ee.Algorithms.String(date.get('year'))
    doy = ee.Algorithms.String(date.getRelative('day', 'year'))
    path = ee.Algorithms.String(img.getNumber('WRS_PATH').toInt())
    row = ee.Algorithms.String(img.getNumber('WRS_ROW').toInt())
    yearDoy = year.cat(doy).cat(path).cat(row)
    return img.set({'imgID': yearDoy,
        'path': path,
        'row': row
    })

def getFootprint(img):
    return ee.Geometry.Polygon(ee.Geometry(img.get('system:footprint')).coordinates())

def filterBounds(aoi):
    return ee.Filter.bounds(aoi)

# def getCfmask(img):
#     cloudShadowBitMask = 1 << 3
#     cloudsBitMask = 1 << 5
#     qa = img.select('pixel_qa')
#     mask = qa.bitwiseAnd(cloudShadowBitMask) \
#         .eq(0) \
#         .And(qa.bitwiseAnd(cloudsBitMask).eq(0))
#     return mask

# def applyCfmask(img):
#     mask = getCfmask(img)
#     return img.updateMask(mask)


def getCfmask(img):
  return img.select('QA_PIXEL').bitwiseAnd(int('11111', 2)).eq(0)

def scaleMask(img):
  def getFactorImg(factorNames):
    factorList = img.toDictionary().select(factorNames).values()
    return ee.Image.constant(factorList)

  scaleImg = getFactorImg(['REFLECTANCE_MULT_BAND_.'])
  offsetImg = getFactorImg(['REFLECTANCE_ADD_BAND_.'])
  scaled = (img.select('SR_B.').multiply(scaleImg).add(offsetImg)
    .multiply(10000).round().int16())

  return (img.addBands(scaled, None, True)
    .select('SR_B.')
    .updateMask(getCfmask(img)))




def viewWrs1Col(params):
    granuleGeom = msslib['getWrs1GranuleGeom'](params['wrs1'])
    params['aoi'] = ee.Geometry(granuleGeom.get('centroid'))
    params['wrs'] = '1'
    mssDnCol = msslib['getCol'](params) \
        .filter(ee.Filter.eq('pr', params['wrs1']))
    msslib['viewThumbnails'](mssDnCol, None)

def getMedoid(col, bands, parallelScale=1):
    col = col.select(bands)
    median = col.reduce(ee.Reducer.median(), parallelScale)

    def mapFun(img):
        dif = ee.Image(img).subtract(median).pow(ee.Image.constant(2))
        return dif.reduce(ee.Reducer.sum()).addBands(img)

    difFromMedian = col.map(mapFun)
    bandNames = difFromMedian.first().bandNames()
    nBands = bandNames.length()
    bandsPos = ee.List.sequence(1, nBands.subtract(1))
    bandNamesSub = bandNames.slice(1)
    return (difFromMedian.reduce(ee.Reducer.min(nBands), parallelScale)
            .select(bandsPos, bandNamesSub))

def getRefImg(params):
    granuleGeoms = msslib['getWrs1GranuleGeom'](params['wrs1'])
    centroid = ee.Geometry(granuleGeoms.get('centroid'))
    bounds = ee.Geometry(granuleGeoms.get('bounds'))

    refCol = msslib['getCol']({
        'aoi': bounds,
        'wrs': '2',
        'yearRange': [1983, 1987], # NOTE: Use five early years, want good coverage, but near to MSS WRS-1 window.
        'doyRange': params['doyRange'],
    }).map(addTmToMssJoinId)

    tmCol = getTmWrs2Col(bounds) \
        .filterDate('1983-01-01', '1988-01-01') \
        .map(addTmToMssJoinId)

    def cloudMask(img):
        mask = getCfmask(ee.Image(img.get('coincidentTmMss')))
        imgToa = msslib['calcToa'](img)
        return imgToa.updateMask(mask)

    coincident = coincidentTmMssCol(refCol, tmCol).map(cloudMask)
    return msslib['addTc'](msslib['addNdvi'](getMedoid(coincident, ['green', 'red', 'red-edge', 'nir']))) \
        .select(['green', 'red', 'red-edge', 'nir', 'ndvi', 'tcb', 'tcg', 'tca']) \
        .set('bounds', bounds)

def exportMssRefImg(params):
    print('Exporting MSS 2nd Gen reference image, please wait.')
    refImg = getRefImg(params)
    task = ee.batch.Export.image.toAsset(**{
        'image': refImg,
        'description': 'MSS-reference-image',
        'assetId': params['baseDir'] + '/ref',
        'region': ee.Geometry(refImg.get('bounds')),
        'scale': 60,
        'crs': params['crs'],
        'maxPixels': 1e13
    })
    task.start()
    return task

# def calcRegression(xImg, yImg, xBand, yBand, aoi, scale):
#     constant = ee.Image(1)
#     xVar = xImg.select(xBand)
#     yVar = yImg.select(yBand)
#     imgRegress = ee.Image.cat(constant, xVar, yVar)

#     linearRegression = imgRegress.reduceRegion(**{
#         'reducer': ee.Reducer.robustLinearRegression(**{
#             'numX': 2,
#             'numY': 1
#         }),
#         'geometry': aoi,
#         'scale': scale,
#         'maxPixels': 1e13
#     })

#     coefList = ee.Array(linearRegression.get('coefficients')).toList()
#     intercept = ee.List(coefList.get(0)).get(0)
#     slope = ee.List(coefList.get(1)).get(0)
#     rmse = ee.Array(linearRegression.get('residuals')).toList().get(0)
#     return ee.Dictionary({'slope': slope, 'intercept': intercept, 'rmse': rmse})

# def applyCoef(img, band, coef):
#     coef = ee.Dictionary(coef)
#     return img.select(band) \
#         .multiply(ee.Image.constant(coef.getNumber('slope'))) \
#         .add(ee.Image.constant(coef.getNumber('intercept')))

# def getSampleImg(img, ref, band):
#     dif = img.select(band) \
#         .subtract(ref.select(band)).rename('dif')

#     difThresh = dif.reduceRegion(**{
#         'reducer': ee.Reducer.percentile(**{
#             'percentiles': [40, 60],
#             'maxRaw': 1000000,
#             'maxBuckets': 1000000,
#             'minBucketWidth': 0.00000000001
#         }),
#         'geometry': img.geometry(),
#         'scale': 60,
#         'maxPixels': 1e13
#     })

#     mask = dif.gt(difThresh.getNumber('dif_p40')) \
#         .And(dif.lt(difThresh.getNumber('dif_p60')))

#     return img.updateMask(mask)

# def correctMssImg(img):
#     ref = ee.Image(img.get('ref_img'))
#     # ref = ee.Image(params['baseDir'] + '/ref')
#     granuleGeoms = msslib['getWrs1GranuleGeom'](img.getString('pr'))
#     granule = ee.Feature(granuleGeoms.get('granule')).geometry()

#     greenCoef = calcRegression(getSampleImg(img, ref, 'green'), ref, 'green', 'green', granule, 60)
#     redCoef = calcRegression(getSampleImg(img, ref, 'red'), ref, 'red', 'red', granule, 60)
#     nirCoef = calcRegression(getSampleImg(img, ref, 'nir'), ref, 'nir', 'nir', granule, 60)
#     ndviCoef = calcRegression(getSampleImg(img, ref, 'ndvi'), ref, 'ndvi', 'ndvi', granule, 60)
#     tcbCoef = calcRegression(getSampleImg(img, ref, 'tcb'), ref, 'tcb', 'tcb', granule, 60)
#     tcgCoef = calcRegression(getSampleImg(img, ref, 'tcg'), ref, 'tcg', 'tcg', granule, 60)
#     tcaCoef = calcRegression(getSampleImg(img, ref, 'tca'), ref, 'tca', 'tca', granule, 60)

#     return ee.Image(ee.Image.cat(
#         applyCoef(img, 'green', greenCoef).toFloat(),
#         applyCoef(img, 'red', redCoef).toFloat(),
#         applyCoef(img, 'nir', nirCoef).toFloat(),
#         applyCoef(img, 'ndvi', nirCoef).toFloat(),
#         applyCoef(img, 'tcb', tcbCoef).toFloat(),
#         applyCoef(img, 'tcg', tcgCoef).toFloat(),
#         applyCoef(img, 'tca', tcaCoef).toFloat()) \
#         .rename(['green', 'red', 'nir', 'ndvi', 'tcb', 'tcg', 'tca']) \
#         .copyProperties(img, img.propertyNames()))

# Makes WRS1 images match the MSS WRS2 reference image using Random Forest regression
# with stratified training sample based on 500 points in 6 categories of TCA defined by
# percentiles
def correctMssImg2(img):
    ref = ee.Image(img.get('ref_img'))
    granuleGeoms = msslib['getWrs1GranuleGeom'](img.getString('pr'))
    granule = ee.Feature(granuleGeoms.get('granule')).geometry()

    img = img.select(ref.bandNames())
    dif = img.subtract(ref).pow(ee.Image.constant(2)).reduce(ee.Reducer.sum())

    difThresh = dif.reduceRegion(**{
        'reducer': ee.Reducer.percentile(**{
            'percentiles': [10],
            'maxRaw': 1000000,
            'maxBuckets': 1000000,
            'minBucketWidth': 0.00000000001
        }),
        'geometry': granule,
        'scale': 60,
        'maxPixels': 1e13
    })

    mask = dif.lt(difThresh.getNumber('sum'));

    tca = img.select('tca')
    tcaGood = tca.gt(0).And(tca.lt(45))
    ndviGood = img.select('ndvi').gt(0)
    tcaMasked = tca.updateMask(mask)#.updateMask(tcaGood).updateMask(ndviGood)
    breaks = [5, 15, 30, 50, 70, 85, 95]
    breakNames = ee.List(breaks).map(lambda num: ee.Number(num).format('%02d'))
    tcaBreaks = tcaMasked.reduceRegion(**{
        'reducer': ee.Reducer.percentile(**{
            'percentiles': breaks,
            'maxRaw': 1000000,
            'maxBuckets': 1000000,
            'minBucketWidth': 0.00000000001,
            'outputNames': breakNames
        }),
        'geometry': granule,
        'scale': 60,
        'maxPixels': 1e13
    })

    breakNames = breakNames.map(lambda i: ee.String('tca_').cat(ee.String(i)))

    for i in range(0, len(breaks)+1):
      if i is 0:
        classImg = tcaMasked.where(
            tcaMasked.lt(tcaBreaks.getNumber(breakNames.get(i))), i)
      elif i is len(breaks):
        classImg = classImg.where(
            tcaMasked.gte(tcaBreaks.getNumber(breakNames.get(i-1))), i)
      else:
        classImg = classImg.where(
            tcaMasked.gte(tcaBreaks.getNumber(breakNames.get(i-1))).And(
            tcaMasked.lt(tcaBreaks.getNumber(breakNames.get(i)))), i)

    sampImg = img.addBands(classImg.rename('class').byte()).addBands(ref)
    sample = sampImg.stratifiedSample(**{
      'numPoints': 500,
      'classBand': 'class',
      'region': granule,
      'scale': 60,
    })

    def predictBand(img, samp, targetBand, inputBands):
      trainedClassifier = ee.Classifier.smileRandomForest(10).train(**{
          'features': samp,
          'classProperty': targetBand+'_1',
          'inputProperties': inputBands
      }).setOutputMode('REGRESSION')
      return img.classify(trainedClassifier).rename(targetBand)

    bandNames = ['green', 'red', 'red-edge', 'nir', 'ndvi', 'tcb', 'tcg', 'tca']
    bands = []
    for band in bandNames:
      bands.append(predictBand(img, sample, band, bandNames))

    return ee.Image(bands).copyProperties(img, img.propertyNames()) # TODO copy only specific properties


# Adds TC and NDVI bands to MSS images and applies QA and MSScvm masks to MSS images
def prepMss(img):
    toa = msslib['calcToa'](img)
    toaAddBands = msslib['addTc'](msslib['addNdvi'](toa))
    toaAddBandsMask = msslib['applyQaMask'](toaAddBands)
    return msslib['applyMsscvm'](toaAddBandsMask)
            #.select(['green', 'red', 'nir', 'ndvi', 'tcb', 'tcg', 'tca'])
            #.multiply(ee.Image([1e4, 1e4, 1e4, 1e4, 1e4, 1e4, 1e2]))
            #.round().toShort().copyProperties(img, img.propertyNames()))

# def processMssWrs1Img(img):
#     toaAddBandsMsscvmMask = prepMss(img)
#     corrected = correctMssImg(toaAddBandsMsscvmMask)
#     return corrected

def processMssWrs1Img2(img):
    toaAddBandsMsscvmMask = prepMss(img)
    corrected = correctMssImg2(toaAddBandsMsscvmMask)
    return corrected

def scaleMssToInt16(img):
  scale = ee.Image([1e4, 1e4, 1e4, 1e4, 1e4, 1e4, 1e4, 1e2])
  return (ee.Image(img.select(['green', 'red', 'red-edge', 'nir', 'ndvi', 'tcb', 'tcg', 'tca'])
          .multiply(scale).round().toShort()
          .copyProperties(img, img.propertyNames()))) # TODO copy only properties needed

# Makes annual MSS image composites based on dates and regions given in params
# The individual MSS WRS1 images are corrected to match the MSS WRS2 reference image
# The images that come out are scaled to int16
# The bands are ['green', 'red', 'red-edge', 'nir', 'ndvi', 'tcb', 'tcg', 'tca']
# Images are 60 meters clipped to the WRS1 tile specified in the params dictionary
# 1983 is MSS WRS2, it is not made to match the reference image
# The export is a single image with bands for each image labeled by year and band name
def processMssWrs1Imgs(params):
    print('Exporting annual MSS composites that match the MSS 2nd Gen reference image, please wait.')
    granuleGeom = msslib['getWrs1GranuleGeom'](params['wrs1'])
    geom = ee.Feature(granuleGeom.get('granule')).geometry()
    params['aoi'] = ee.Geometry(granuleGeom.get('centroid'))
    params['wrs'] = '1'

    def setRefImg(img):
        return img.set('ref_img', ee.Image(params['baseDir'] + '/ref'))

    mssCol = (msslib['getCol'](params)
        .filter(ee.Filter.eq('pr', params['wrs1'])))
        
    mss1983 = msslib['getCol']({
        'aoi': geom,
        'wrs': '2',
        'yearRange': [1983, 1983],
        'doyRange': params['doyRange']
    })

    mssCol = mssCol.merge(mss1983).map(setRefImg)

    dummy = (ee.Image([0, 0, 0, 0, 0, 0, 0, 0]).selfMask().toShort()
        .rename(['green', 'red', 'red-edge', 'nir', 'ndvi', 'tcb', 'tcg', 'tca']))
    
    imgs = []
    for y in range(1972, 1984):
        print('Exporting year:', y)
        yrCol = mssCol.filter(ee.Filter.eq('year', y))
        n_imgs = yrCol.size().getInfo()
        if (n_imgs == 0):
            print('  no images, exporting placeholder')
            yearImg = dummy
        else:
            if y != 1983:
                yrCol = yrCol.map(processMssWrs1Img2)
                parallelScale = 1
            else:
                yrCol = yrCol.map(prepMss)  # NOTE: not 1983 normalized to ref image
                parallelScale = 4
            yearImg = getMedoid(yrCol, ['green', 'red', 'red-edge', 'nir', 'ndvi', 'tcb', 'tcg', 'tca'], parallelScale)
        
        yearImg = scaleMssToInt16(yearImg).set('year', y)
        yearImg = appendYearToBandnames(yearImg)
        imgs.append(yearImg)

    outImg = appendIdToBandnames(ee.ImageCollection(imgs).toBands())
    task = ee.batch.Export.image.toAsset(**{
        'image': outImg.clip(geom),
        'description': 'MSS_WRS1_to_WRS2_stack',
        'assetId': params['baseDir'] + '/MSS_WRS1_to_WRS2_stack',
        'region': geom,
        'scale': 60,
        'crs': params['crs']
    })
    task.start()
    return task



# def correctMss1983(params):
#     aoi = ee.Feature(
#         msslib['getWrs1GranuleGeom'](params['wrs1']).get('granule')).geometry()
#     mssCol = msslib['getCol']({
#         'aoi': aoi,
#         'wrs': '2',
#         'yearRange': [1983, 1983],
#         'doyRange': params['doyRange']
#     }).map(prepMss)
    
#     mssCol1983 = (mssCol.map(correctMssImg_doit)
#       .map(lambda img: img.resample('bicubic')))
    
#     outImg = getMedoid(mssCol1983, ['blue', 'green', 'red', 'nir', 'ndvi', 'tcb', 'tcg', 'tcw', 'tca']) \
#         .round().toShort().clip(aoi).set({
#             'dummy': False,
#             'year': 1983,
#             'system:time_start': ee.Date.fromYMD(1983, 1, 1)
#         })

#     task = ee.batch.Export.image.toAsset(**{
#         'image': outImg,
#         'description': 'WRS1_to_TM_1983',
#         'assetId': params['baseDir'] + '/WRS1_to_TM/' + '1983',
#         'region': aoi,
#         'scale': 30,
#         'crs': params['crs'],
#         'maxPixels': 1e13
#     })
#     task.start()
#     return [task]
#     # imgs = mssColToTm.aggregate_array('system:index').getInfo()
#     # tasks = []
#     # for i in range(0, len(imgs)):
#     #   fname = '1983_' + str(i).zfill(2)
#     #   print(fname)
#     #   thisImg = mssColToTm.filter(ee.Filter.eq('system:index', imgs[i])).first()

#     #   task = ee.batch.Export.image.toAsset(**{
#     #     'image': thisImg.clip(aoi),
#     #     'description': fname,
#     #     'assetId': params['baseDir'] + '/mss_1983_col/' + fname,
#     #     'region': aoi,
#     #     'scale': 60,  # should this be 30
#     #     'crs': params['crs'],
#     #     'maxPixels': 1e13
#     #   })
#     #   task.start()
#     #   tasks.append(task)
    
#     # return tasks

def renameOli(img):
    return img.select(
        ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'],  # , 'QA_PIXEL'
        ['blue', 'green', 'red', 'nir', 'swir1', 'swir2']) # , 'pixel_qa']

def renameTm(img):
    return img.select(
        ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'], # , 'QA_PIXEL'
        ['blue', 'green', 'red', 'nir', 'swir1', 'swir2']) # , 'pixel_qa'

def tmAddIndices(img):
    b = ee.Image(img).select(['blue', 'green', 'red', 'nir', 'swir1', 'swir2'])
    brt_coeffs = ee.Image.constant([0.2043, 0.4158, 0.5524, 0.5741, 0.3124, 0.2303])
    grn_coeffs = ee.Image.constant([-0.1603, -0.2819, -0.4934, 0.7940, -0.0002, -0.1446])
    wet_coeffs = ee.Image.constant([0.0315, 0.2021, 0.3102, 0.1594, -0.6806, -0.6109])
    brightness = b.multiply(brt_coeffs).reduce(ee.Reducer.sum()).round().toShort().rename('tcb')
    greenness = b.multiply(grn_coeffs).reduce(ee.Reducer.sum()).round().toShort().rename('tcg')
    wetness = b.multiply(wet_coeffs).reduce(ee.Reducer.sum()).round().toShort().rename('tcw')
    angle = (greenness.divide(brightness)).atan().multiply(180 / math.pi).multiply(100).round().toShort().rename('tca')
    ndvi = img.normalizedDifference(['nir', 'red']).multiply(1000).round().toShort().rename('ndvi')
    tc = ee.Image.cat(ndvi, brightness, greenness, wetness, angle)#.rename(['ndvi', 'tcb', 'tcg', 'tcw' 'tca'])
    return img.addBands(tc)

def gatherTmCol(params):
    granuleGeom = msslib['getWrs1GranuleGeom'](params['wrs1'])
    aoi = ee.Feature(granuleGeom.get('granule')).geometry()
    dateFilter = ee.Filter.calendarRange(params['doyRange'][0], params['doyRange'][1], 'day_of_year')
    startDate = ee.Date.fromYMD(params['yearRange'][0], 1, 1)
    endDate = startDate.advance(1, 'year')
    oli2Col = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2') \
        .filterBounds(aoi).filterDate(startDate, endDate).filter(dateFilter).map(prepOli)
    oliCol = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2') \
        .filterBounds(aoi).filterDate(startDate, endDate).filter(dateFilter).map(prepOli)
    etmCol = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2') \
        .filterBounds(aoi).filterDate(startDate, endDate).filter(dateFilter).map(prepTm)
    tm5Col = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2') \
        .filterBounds(aoi).filterDate(startDate, endDate).filter(dateFilter).map(prepTm)
    tm4Col = ee.ImageCollection('LANDSAT/LT04/C02/T1_L2') \
        .filterBounds(aoi).filterDate(startDate, endDate).filter(dateFilter).map(prepTm)
    return ee.ImageCollection(ee.FeatureCollection([tm4Col, tm5Col, etmCol, oliCol, oli2Col]).flatten())

def prepOli(img):
    orig = img
    img = scaleMask(img)
    img = renameOli(img)
    img = tmAddIndices(img)
    return ee.Image(img.copyProperties(orig, orig.propertyNames())) # TODO: only copy the needed properties

def prepTm(img):
    orig = img
    img = scaleMask(img)
    img = renameTm(img)
    img = tmAddIndices(img)
    return ee.Image(img.copyProperties(orig, orig.propertyNames())) # TODO: only copy the needed properties

def getCoincidentTmMssCol(params):
    aoi = ee.Feature(
    msslib['getWrs1GranuleGeom'](params['wrs1']).get('granule')).geometry()
    mssCol = msslib['getCol']({
        'aoi': aoi,
        'wrs': '2',
        'doyRange': params['doyRange'],
        'excludeIds': params['excludeIds']
    }) \
    .map(addTmToMssJoinId)

    tmCol = getTmWrs2Col(aoi).map(addTmToMssJoinId)
    coincident = coincidentTmMssCol(mssCol, tmCol)
    return coincident

# Gets a sample of pixels from coincident MSS and TM images. The TM image ID
# is a property of the input MSS image. The bands of TM image are added to the
# MSS image and then sampled using a stratified class band based on MSS TCA
# percentile bins.
def getMsstoTmStratSamp(img):
    xImg = scaleMssToInt16(msslib['addTc'](msslib['addNdvi'](msslib['calcToa'](img))))
    yImg = prepTm(ee.Image(xImg.get('coincidentTmMss')))
    granule = ee.Feature(ee.FeatureCollection('users/jstnbraaten/wrs/wrs2_descending_land') \
        .filter(ee.Filter.eq('PR', xImg.getString('pr'))).first()).geometry()

    tca = xImg.select('tca')
    #tcaGood = tca.gt(0).And(tca.lt(45))
    ndviGood = xImg.select('ndvi').gt(-500)
    tcaMasked = tca#.updateMask(tcaGood).updateMask(ndviGood)
    breaks = [5, 15, 30, 50, 70, 85, 95]
    breakNames = ee.List(breaks).map(lambda num: ee.Number(num).format('%02d'))
    tcaBreaks = tcaMasked.reduceRegion(**{
        'reducer': ee.Reducer.percentile(**{
            'percentiles': breaks,
            'maxRaw': 1000000,
            'maxBuckets': 1000000,
            'minBucketWidth': 0.00000000001,
            'outputNames': breakNames
        }),
        'geometry': granule,
        'scale': 60,
        'maxPixels': 1e13
    })

    breakNames = breakNames.map(lambda i: ee.String('tca_').cat(ee.String(i)))

    for i in range(0, len(breaks)+1):
      if i is 0:
        classImg = tcaMasked.where(
            tcaMasked.lt(tcaBreaks.getNumber(breakNames.get(i))), i)
      elif i is len(breaks):
        classImg = classImg.where(
            tcaMasked.gte(tcaBreaks.getNumber(breakNames.get(i-1))), i)
      else:
        classImg = classImg.where(
            tcaMasked.gte(tcaBreaks.getNumber(breakNames.get(i-1))).And(
            tcaMasked.lt(tcaBreaks.getNumber(breakNames.get(i)))), i)

    sampImg = xImg.addBands(classImg.rename('class').byte()).addBands(yImg)
    sample = sampImg.stratifiedSample(**{
      'numPoints': 25,  # Is this value enough, could it be higher, check the output table
      'classBand': 'class',
      'region': granule,
      'scale': 60,
      'geometries': True
    })   

    return (sample.map(lambda f: f.set('constant', 1))
      .copyProperties(img, ['imgID', 'year', 'path', 'row', 'pr'])) # TODO: these properties are not present in output, probably not in input img

# Creates a sample of pixels from coincident TM and MSS images to use a training
# sample in a random forest regression classifier. The sample is stratified on MSS
# TCA from 8(?) percentile bins. 
def exportMss2TmCoefCol(params):
    print('Exporting MSS-to-TM model training sample, please wait.')
    col = getCoincidentTmMssCol(params)
    sample = col.map(getMsstoTmStratSamp).flatten()
    task = ee.batch.Export.table.toAsset(**{
        'collection': sample,
        'description': 'mss_to_tm_coef_fc',
        'assetId': params['baseDir'] + '/mss_to_tm_coef_fc'
    })
    print('Exporting mss_to_tm_coef_fc')
    task.start()
    return task

def predictBand_doit(sample, img, targetBand, outName):
  bands = ['green', 'red', 'red-edge', 'nir', 'ndvi', 'tcb', 'tcg', 'tca']
  trainedClassifier = ee.Classifier.smileRandomForest(10).train(**{
    'features': sample,
    'classProperty': targetBand,
    'inputProperties': bands
  }).setOutputMode('REGRESSION')
  return img.classify(trainedClassifier).rename(outName).round().toShort() 

def correctMssImg_doit(img):
  print('1055')
  print('params', params)
  print("params['baseDir'] + '/mss_to_tm_coef_fc'", params['baseDir'] + '/mss_to_tm_coef_fc')
  sample = ee.FeatureCollection(params['baseDir'] + '/mss_to_tm_coef_fc')
  print('1059')
  targetBands = ['blue', 'green_1', 'red_1', 'nir_1', 'swir1', 'swir2', 'ndvi_1', 'tcb_1', 'tcg_1', 'tcw', 'tca_1'] #['blue', 'green_1', 'red_1', 'nir_1', 'ndvi_1', 'tcb_1', 'tcg_1', 'tca_1']
  outBands = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'ndvi', 'tcb', 'tcg', 'tcw', 'tca'] # ['blue', 'green', 'red', 'nir', 'ndvi', 'tcb', 'tcg', 'tca']
  # bands = []
  bands = [None for b in range(0, len(outBands))]
  for i in range(0, len(outBands)):
    # print(i)
    # print(targetBands[i])
    print(outBands[i])
    band = predictBand_doit(sample, img, targetBands[i], outBands[i])
    # bands.append(band)
    bands[i] = band

  return ee.Image(bands).copyProperties(img, img.propertyNames())


def mssStackToCol(mssStackPath):
    imgStack = ee.Image(mssStackPath)
    imgList = []
    for y in range(1972, 1984):
        img = getImgYearFromStack(imgStack, y)
        imgList.append(img)
    return ee.ImageCollection(imgList)

def getFinalCorrectedMssCol(imgStackPath):
    imgStack = ee.Image(imgStackPath)
    imgList = []
    for y in range(1972, 1984):
        img = getImgYearFromStack(imgStack, y)
        imgList.append(img)

    return appendIdToBandnames(ee.ImageCollection(imgList)
                               .map(correctMssImg_doit)
                               .map(appendYearToBandnames)
                               .toBands())


def exportFinalCorrectedMssCol(params):
    print('Exporting annual MSS composites that match TM, please wait.')
    mssCol = mssStackToCol(params['baseDir'] + '/MSS_WRS1_to_WRS2_stack')
    outImg = appendIdToBandnames(mssCol
                                 .map(correctMssImg_doit)
                                 .map(appendYearToBandnames)
                                 .toBands())

    granuleGeom = msslib['getWrs1GranuleGeom'](params['wrs1'])
    geom = ee.Feature(granuleGeom.get('granule')).geometry()
    
    print('Exporting MSS to TM stack:')
    task = ee.batch.Export.image.toAsset(**{
        'image': outImg.resample('bicubic').clip(geom),
        'description': 'WRS1_to_TM_stack',
        'assetId': params['baseDir'] + '/WRS1_to_TM_stack',
        'region': geom,
        'scale': 30,
        'crs': params['crs'],
        'maxPixels': 1e13
    })
    task.start()
    return task


def exportMss1983(params):
    aoi = ee.Feature(
        msslib['getWrs1GranuleGeom'](params['wrs1']).get('granule')).geometry()
    mssCol1983 = ee.ImageCollection(params['baseDir'] + '/mss_1983_col') \
      .map(lambda img: img.resample('bicubic'))
    
    outImg = (getMedoid(mssCol1983, ['blue', 'green', 'red', 'nir', 'swir1', 'ndvi', 'tcb', 'tcg', 'tca'])  #['blue', 'green', 'red', 'nir', 'ndvi', 'tcb', 'tcg', 'tca']) \
        .set({
            'dummy': False,
            'year': 1983,
            'system:time_start': ee.Date.fromYMD(1983, 1, 1)
        })).toShort().clip(aoi)

    task = ee.batch.Export.image.toAsset(**{
        'image': outImg,
        'description': 'WRS1_to_TM' + '1983',
        'assetId': params['baseDir'] + '/WRS1_to_TM/' + '1983',
        'region': aoi,
        'scale': 30,
        'crs': params['crs'],
        'maxPixels': 1e13
    })
    task.start()
    return [task]


def runLt(params):  #exportTmComposites(params):
    def add_systime(img):
      millis = ee.Date.fromYMD(img.getNumber('year'), 6 ,1).millis()
      date = ee.Date(millis).format('YYYY-MM-dd')
      return img.set({'system:time_start': millis, 'date': date})

    granuleGeom = msslib['getWrs1GranuleGeom'](params['wrs1'])
    geom = ee.Feature(granuleGeom.get('granule')).geometry()

    mssCol = mssStackToCol(params['baseDir'] + '/WRS1_to_TM_stack').map(add_systime) # TODO: these images do not have a system:time_start - they should have one - should not need to add
    # for y in range(1972, 2022):
    #   print(mssCol.filter(ee.Filter.eq('year', y)).first().bandNames().getInfo())

    dummyImg = mssCol.first()
    dummyImg = dummyImg.updateMask(dummyImg.mask().multiply(0))

    todaysDate = date.today()
    # tmList = []

    def getTmYearImg(year):
        params_ = copy.deepcopy(params)
        params_['yearRange'] = [year, year]
        bands = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'ndvi', 'tcb', 'tcg', 'tcw', 'tca']  # ['blue', 'green', 'red', 'nir', 'ndvi', 'tcb', 'tcg', 'tca'])
        parallelScale = 8
        tmCol = gatherTmCol(params_)
        date = ee.Date(ee.Date.fromYMD(year, 6 ,1).millis())        
        return ee.Image(ee.Algorithms.If(tmCol.size(), getMedoid(tmCol, bands, parallelScale), dummyImg)).set({
            'system:time_start': date.millis(),
            'year': date.get('year'), 
            'date': date.format('YYYY-MM-dd')
            })
    tmYearRange = ee.List.sequence(1984, 2022) # NOTE: end year is set here. -todaysDate.year
    tmYearImgList = tmYearRange.map(getTmYearImg)
    tmYearImgCol = ee.ImageCollection.fromImages(tmYearImgList)

        



    # for y in range(1984, 2022):  # todaysDate.year + 1
    #     params['yearRange'] = [y, y]
    #     thisYearImg = (getMedoid(
    #         gatherTmCol(params), ['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'ndvi', 'tcb', 'tcg', 'tcw', 'tca'])  # ['blue', 'green', 'red', 'nir', 'ndvi', 'tcb', 'tcg', 'tca'])
    #         .set({
    #             'system:time_start': ee.Date.fromYMD(y, 1 ,1).millis(),
    #             'year': y
    #         }))
    #     tmList.append(thisYearImg)


    
    def addSegBand(img):
      outImg = img.select('ndvi').multiply(-1).rename('seg')
      return outImg.addBands(img).copyProperties(img, ['system:time_start', 'year', 'date'])

    #tmCol = ee.ImageCollection(tmList)
    #col = mssCol.merge(tmCol).map(addSegBand).sort('year')

    col = mssCol.merge(tmYearImgCol).map(addSegBand).sort('year')
    years = col.aggregate_array('year')#.getInfo()
    # millis = col.aggregate_array('system:time_start').getInfo()

    # print('all years', years)
    # print('all millis', millis)



    # for y in range(1972, 2022):
    #   print(col.filter(ee.Filter.eq('year', y)).first().bandNames().getInfo())

    # Deal with missing years (check is other processes are filling bads years, I think MSS is)
    # year_ends = col.aggregate_array('year').reduce(ee.Reducer.minMax()).getInfo()
    # print('year_ends', year_ends)
    # Make a dummy image for missing years. 
    # bandNames = ee.List(['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2', 'pixel_qa'])
    # fillerValues = ee.List.repeat(0, bandNames.size())
    # dummyImg = ee.Image.constant(fillerValues).rename(bandNames).selfMask().int16()


    # dummies = []
    # years = list(range(year_ends['min'], year_ends['max']+1))
    # print('years', years)
    # for y in range(year_ends['min'], year_ends['max']+1):
    #   print('Checking year', y)
    #   colSize = col.filter(ee.Filter.eq('year', y)).size().getInfo()
    #   if colSize == 0:
    #     print(y, ' is a dummy')
    #     dummies.append(dummyImg.set({
    #         'system:time_start': ee.Date.fromYMD(y, 6 ,1).millis(),
    #         'year': y
    #         }))

    # col = col.merge(ee.ImageCollection(dummies)).sort('year')



    # col = col.map(add_systime)
    # print('dates', col.aggregate_array('date').getInfo())
    lt = ee.Algorithms.TemporalSegmentation.LandTrendr(**{
        'timeSeries': col.select(['seg'] + params['ltParams']['ftvBands']),
        'maxSegments': params['ltParams']['maxSegments'],
        'spikeThreshold': params['ltParams']['spikeThreshold'],
        'vertexCountOvershoot': params['ltParams']['vertexCountOvershoot'],
        'preventOneYearRecovery': params['ltParams']['preventOneYearRecovery'],
        'recoveryThreshold': params['ltParams']['recoveryThreshold'],
        'pvalThreshold': params['ltParams']['pvalThreshold'],
        'bestModelProportion': params['ltParams']['bestModelProportion'],
        'minObservationsNeeded': params['ltParams']['minObservationsNeeded']
    }).set({'years': years})

    return lt.int16()




def exportLt(params):
    print('Exporting LandTrendr segmentation and FTV image array, please wait.')
    lt = runLt(params)
    years = ee.Image(lt).get('years').getInfo()
    yearsStr = ['yr_' + str(year) for year in years]
    granuleGeom = msslib['getWrs1GranuleGeom'](params['wrs1'])
    geom = ee.Feature(granuleGeom.get('granule')).geometry()

    # tasks = []
    # for band in params['ltParams']['ftvBands']:
    #   print(getPaths(params['projectDir'], params['wrs1'])['fit_collection']+ '/' + band + '_fit')
    #   ftv_img = lt.select([band + '_fit']).arrayFlatten([yearsStr]).toShort().set('band', band)
    #   task = ee.batch.Export.image.toAsset(**{
    #       'image': ftv_img.clip(geom),
    #       'description': band + '_fit',
    #       'assetId': getPaths(params['projectDir'], params['wrs1'])['fit_collection'] + '/' + band + '_fit',
    #       'region': geom,
    #       'scale': 30,
    #       'crs': params['crs'],
    #       'maxPixels': 1e13
    #   })
    #   task.start()
    #   tasks.append(task)
    otherBands = ['rmse'] + [f'{band}_fit' for band in params['ltParams']['ftvBands']]
    ltlt = (lt.select('LandTrendr').arrayPad([4, len(years)], 0)
              .addBands(lt.select(otherBands)).int16())

    ltltBands = ['LandTrendr'] + otherBands
    pyramidingPolicy = {}
    for band in ltltBands:
        if band in ['LandTrendr', 'rmse']:
            pyramidingPolicy[band] = 'sample'
        else:
            pyramidingPolicy[band] = 'mean'

    print(params['baseDir'] + '/landtrendr')
    task = ee.batch.Export.image.toAsset(**{
        'image': ltlt.clip(geom),
        'pyramidingPolicy': pyramidingPolicy,
        'description': 'landtrendr',
        'assetId': params['baseDir'] + '/landtrendr',
        'region': geom,
        'scale': params['ltParams']['scale'],
        'crs': params['crs'],
        'maxPixels': 1e13,
        'shardSize': 32
    })
    task.start()
    # tasks.append(task)

    return task






def appendYearToBandnames(img):
  img = ee.Image(img)
  year = img.getNumber('year').format()
  names = img.bandNames()
  namesYear = names.map(lambda name: ee.String(year).cat(ee.String('_')).cat(ee.String(name)))
  return img.rename(namesYear)

def appendIdToBandnames(img):
  def processName(band):
    parts = ee.String(band).split('_')
    id = ee.String('img').cat(ee.Number.parse(parts.getString(0)).format('%02d'))
    return parts.splice(0, 1, [id]).join('_')
  
  img = ee.Image(img)
  names = img.bandNames()
  namesId = names.map(processName)
  return img.rename(namesId)

def getImgYearFromStack(img, year):
  search = 'img.._' + str(year) + '.*'
  theseBands = img.select(search)
  namesLong = theseBands.bandNames()
  namesShort = namesLong.map(lambda name: ee.String(name).split('_').get(-1))
  return theseBands.rename(namesShort).set('year', year)


# def getColForLandTrendrOnTheFly(params):
#     mssCol = getFinalCorrectedMssCol(params)

#     tmCol = ee.ImageCollection([])
#     for y in range(1983, 2022):  # make this the current year
#         params['yearRange'] = [y, y]
#         thisYearCol = getMedoid(gatherTmCol(params), ['blue', 'green', 'red', 'nir', 'ndvi', 'tcb', 'tcg', 'tca']) \
#             .set('system:time_start', ee.Date.fromYMD(y, 1 ,1).millis())
#         tmCol = tmCol.merge(ee.ImageCollection(thisYearCol.toShort()))

#     def prepForLt(img):
#         return img.select('ndvi').multiply(-1).rename('LTndvi') \
#             .addBands(img.select(['blue', 'green', 'red', 'nir', 'ndvi', 'tcb', 'tcg', 'tca'])) \
#             .set('system:time_start', img.get('system:time_start'))

#     combinedCol = mssCol.merge(tmCol).map(prepForLt).sort('system:time_start')

#     return combinedCol



# def getColForLandTrendrFromAsset(params):  # Relies on WRS1_to_TM assets:
#     mssCol = ee.ImageCollection([])
#     for y in range(1972, 1983):
#         img = ee.Image(params['baseDir'] + '/WRS1_to_TM/' + str(y)) \
#             .set('system:time_start', ee.Date.fromYMD(y, 1 ,1).millis())
#         mssCol = mssCol.merge(ee.ImageCollection(img))

#     mss1983 = ee.ImageCollection(
#         getMedoid(correctMssWrs2(params), ['blue', 'green', 'red', 'nir', 'ndvi', 'tcb', 'tcg', 'tca']) \
#         .set('system:time_start', ee.Date.fromYMD(1983, 1 ,1).millis()))

#     tmCol = ee.ImageCollection([])
#     for y in range(1984, 2020, 1):
#         params['yearRange'] = [y, y]
#         thisYearCol = getMedoid(
#             gatherTmCol(params), ['blue', 'green', 'red', 'nir', 'ndvi', 'tcb', 'tcg', 'tca']) \
#             .set('system:time_start', ee.Date.fromYMD(y, 1 ,1).millis())
#         tmCol = tmCol.merge(ee.ImageCollection(thisYearCol.toShort()))

#     def prepForLt(img):
#         return img.select('ndvi').multiply(-1).rename('LTndvi') \
#             .addBands(img.select(['blue', 'green', 'red', 'nir', 'ndvi', 'tcb', 'tcg', 'tca'])) \
#             .set('system:time_start', img.get('system:time_start'))

#     return mssCol.merge(mss1983).merge(tmCol).map(prepForLt).sort('system:time_start')

# def runLandTrendrMss2Tm(params):
#     ltCol = getColForLandTrendrFromAsset(params)
#     lt = ee.Algorithms.TemporalSegmentation.LandTrendr(**{
#         'timeSeries': ltCol,
#         'maxSegments': 10,
#         'spikeThreshold': 0.7,
#         'vertexCountOvershoot': 3,
#         'preventOneYearRecovery': True,
#         'recoveryThreshold': 0.5,
#         'pvalThreshold': 0.05,
#         'bestModelProportion': 0.75,
#         'minObservationsNeeded': 6
#     })
#     return lt


def checkStatus(taskList):
  status = []
  for i in range(0, len(taskList)):
    status.append(taskList[i].status()['state'] == 'COMPLETED')
  return all(status)


def printStatus(taskList):
  for i in range(0, len(taskList)):
    print(i)
    print(taskList[i].status()['state'])

def monitorTaskStatus(task):
    keep_going = True
    while keep_going:
        time.sleep(60)
        state = task.status()['state']
        if state not in ['UNSUBMITTED', 'READY', 'RUNNING']:
            keep_going = False
    else:
        print(state)

# Code to create and delete assets.
def getPaths(project_dir, wrs_1_granule):
  base_dir = os.path.join(project_dir, wrs_1_granule)
  #mss_wrs1_to_wrs2 = os.path.join(base_dir, 'WRS1_to_WRS2')
  #mss_wrs1_to_tm = os.path.join(base_dir, 'WRS1_to_TM')
  #post_mss = os.path.join(base_dir, 'post_mss')
  #collection = os.path.join(base_dir, 'collection')
  mss_1983_col = os.path.join(base_dir, 'mss_1983_col')
  fit_collection = os.path.join(base_dir, 'fit_collection')
  return {
      'wrs1_scene': base_dir,
      #'mss_wrs1_to_wrs2': mss_wrs1_to_wrs2,
      #'mss_wrs1_to_tm': mss_wrs1_to_tm,
      #'post_mss': post_mss,
      #'collection': collection,
      #'mss_1983_col': mss_1983_col,
      'fit_collection': fit_collection
  }

def createProjectDir(project_dir, wrs_1_granule):
  paths = getPaths(project_dir, wrs_1_granule)
  os.system('earthengine create folder ' + paths['wrs1_scene'])
  #os.system('earthengine create folder ' + paths['mss_wrs1_to_wrs2'])
  #os.system('earthengine create folder ' + paths['mss_wrs1_to_tm'])
  #os.system('earthengine create folder ' + paths['post_mss'])
  #os.system('earthengine create collection ' + paths['collection'])
  #os.system('earthengine create collection ' + paths['mss_1983_col'])
  os.system('earthengine create collection ' + paths['fit_collection'])

# def rmMssWrs1ToWrs2(params):
#   paths = getPaths(params['projectDir'], params['wrs1'])
#   assets = ee.data.listAssets(
#       params={'parent': paths['mss_wrs1_to_wrs2']})['assets']
#   for asset in assets:
#     os.system(' '.join(['earthengine rm', asset['name']]))

# def rmMssWrs1ToTm(params):
#   paths = getPaths(params['projectDir'], params['wrs1'])
#   assets = ee.data.listAssets(
#       params={'parent': paths['mss_wrs1_to_tm']})['assets']
#   for asset in assets:
#     os.system(' '.join(['earthengine rm', asset['name']]))

# def rmPostMss(params):
#   paths = getPaths(params['projectDir'], params['wrs1'])
#   assets = ee.data.listAssets(
#       params={'parent': paths['post_mss']})['assets']
#   for asset in assets:
#     os.system(' '.join(['earthengine rm', asset['name']]))

# def rmCollection(params):
#   paths = getPaths(params['projectDir'], params['wrs1'])
#   assets = ee.data.listImages(
#       params={'parent': paths['collection']})['images']
#   for asset in assets:
#     os.system(' '.join(['earthengine rm', asset['name']]))

def rmFitCollection(params):
  paths = getPaths(params['projectDir'], params['wrs1'])
  assets = ee.data.listImages(
      params={'parent': paths['fit_collection']})['images']
  for asset in assets:
    os.system(' '.join(['earthengine rm', asset['name']]))

def rmMss2TmInfo(params):
  os.system(' '.join(['earthengine rm', os.path.join(params['baseDir'], 'mss_to_tm_coef_fc')]))
#  os.system(' '.join(['earthengine rm', os.path.join(params['baseDir'], 'mss_offset')]))

# def images2Col_doit(path, outCol):
#   assets = ee.data.listAssets(params={'parent': path})['assets']
#   for asset in assets:
#     old = asset['name']
#     year = os.path.basename(old)
#     print('Moving:', year)
#     new = os.path.join(outCol, year)
#     os.system(' '.join(['earthengine rm', new]))
#     os.system(' '.join(['earthengine cp', old, new]))

# def images2Col(params):
#   paths = getPaths(params['projectDir'], params['wrs1'])
#   _images2Col(paths['mss_wrs1_to_tm'], paths['collection'])
#   #_images2Col(paths['post_mss'], paths['collection'])
