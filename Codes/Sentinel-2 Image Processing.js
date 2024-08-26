/*
* Author: Hui Li
* Institute: CSISS GMU
* Date: August. 2024
* e-mail: hli47@gmu.edu
* This program conducts Sentinel2 image collection and composite image calculation.
*/
//method of Sentinel-2 collection. Reference: https://developers.google.com/earth-engine/tutorials/community/sentinel-2-s2cloudless
var S2_cutCldSlw = function(start_date,end_date,boundary,CLOUD_FILTER,CLD_PRB_THRESH,NIR_DRK_THRESH,CLD_PRJ_DIST,BUFFER){
  // Function to add NDVI, time, and constant variables to Sentinel imagery.
  var addVariables = function(image) {
    // Compute time in fractional years since the epoch.
    var date = ee.String(image.get('system:index'));
    // var days = ee.Date(date.slice(0,4).cat('-').cat(date.slice(4,6)).cat('-').cat(date.slice(6,8))).format('DDD');
    var year = date.slice(0,4);
    var month = date.slice(4,6);
    var dateOfMonth = date.slice(6,8);
    //generating day number in the imgcollection
    var days = ee.Number.parse(ee.Date(year.cat('-').cat(month).cat('-').cat(dateOfMonth)).format('DDD'))
                        .add((ee.Number.parse(year).subtract(ee.Number.parse(startYear))).multiply(365))
                        .subtract(initialDays)

    // Return the image with the added bands.
    return image
    // Add an NDVI band.
    .addBands(image.normalizedDifference(['B8', 'B4']).float().rename('NDVI'))
    .set('system:day_start',ee.Number.parse(days));
  };
  /************************************************************/
  // Remove clouds, add variables and filter to the area of interest.
  // Cloud components
  var add_cloud_bands = function(img){
      //Get s2cloudless image, subset the probability band.
      var cld_prb = ee.Image(img.get('s2cloudless')).select('probability');
      //Condition s2cloudless by the probability threshold value.
      var is_cloud = cld_prb.gt(CLD_PRB_THRESH).rename('clouds');
      //Add the cloud probability layer and cloud mask as image bands.
      return img.addBands(ee.Image([cld_prb, is_cloud]));
  }
  //Cloud shadow components
  var add_shadow_bands = function(img){
      //Identify water pixels from the SCL band.
      //var not_water = img.select('SCL').neq(6)
      //Identify dark NIR pixels that are not water (potential cloud shadow pixels).
      var SR_BAND_SCALE = 1e4
      var dark_pixels = img.select('B8').lt(NIR_DRK_THRESH*SR_BAND_SCALE).rename('dark_pixels')//.multiply(not_water)
      //Determine the direction to project cloud shadow from clouds (assumes UTM projection).
      var shadow_azimuth = ee.Number(90).subtract(ee.Number(img.get('MEAN_SOLAR_AZIMUTH_ANGLE')));
      //Project shadows from clouds for the distance specified by the CLD_PRJ_DIST input.
      var cld_proj = (img.select('clouds').directionalDistanceTransform(shadow_azimuth, CLD_PRJ_DIST*10)
          // .reproject({'crs': img.select(0).projection(), 'scale': 100})
          .select('distance')
          .mask()
          .rename('cloud_transform'))
  
      //Identify the intersection of dark pixels with cloud shadow projection.
      var shadows = cld_proj.multiply(dark_pixels).rename('shadows')
      //Add dark pixels, cloud projection, and identified shadows as image bands.
      return img.addBands(ee.Image([dark_pixels, cld_proj, shadows]))
  }
  //Final cloud-shadow mask
  var add_cld_shdw_mask = function(img){
      //Add cloud component bands.
      var img_cloud = add_cloud_bands(img)
      //Add cloud shadow component bands.
      var img_cloud_shadow = add_shadow_bands(img_cloud)
      //Combine cloud and shadow mask, set cloud and shadow as value 1, else 0.
      var is_cld_shdw = img_cloud_shadow.select('clouds').add(img_cloud_shadow.select('shadows')).gt(0)
      //Remove small cloud-shadow patches and dilate remaining pixels by BUFFER input.
      //20 m scale is for speed, and assumes clouds don't require 10 m precision.
      is_cld_shdw = (is_cld_shdw.focalMin(2).focalMax(BUFFER*2/20)
          // .reproject({'crs': img.select([0]).projection(), 'scale': 20})
          .rename('cloudmask'))
      //Add the final cloud-shadow mask to the image.
      return img_cloud_shadow.addBands(is_cld_shdw)
  }
  
  var apply_cld_shdw_mask = function(img){
      //Subset the cloudmask band and invert it so clouds/shadow are 0, else 1.
      var not_cld_shdw = img.select('cloudmask').not()
      //Subset reflectance bands and update their masks, return the result.
      return img.select('B.*').updateMask(not_cld_shdw)
  }
  //mask the water pixel
  var apply_scl_water_mask = function(img){
    var scl = img.select('SCL');
    var wantedPixels = scl.neq(6);
    var targetPixels = scl.eq(4).or(scl.eq(5));
    return img.updateMask(wantedPixels);
  }
  // Import and filter S2 SR.
  var s2_sr_col = ee.ImageCollection('COPERNICUS/S2_HARMONIZED')//S2_SR_HARMONIZED
          .filterDate(start_date, end_date)
          .filterBounds(boundary)
          .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', CLOUD_FILTER))
          .filter(ee.Filter.eq('GENERAL_QUALITY','PASSED'));
  // Import and filter s2cloudless.
  var s2_cloudless_col = ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')
          .filterBounds(boundary)
          .filterDate(start_date, end_date)    
  var s2_sr_cld_col_eval = ee.ImageCollection(ee.Join.saveFirst('s2cloudless').apply({
          'primary': s2_sr_col,
          'secondary': s2_cloudless_col,
          'condition': ee.Filter.equals({
              'leftField': 'system:index',
              'rightField': 'system:index'
          })
      }))
  //finding the initial days of the year    
  var initialDate = ee.String(s2_sr_col.first().get('system:index'));
  var startYear = initialDate.slice(0,4);
  var month = initialDate.slice(4,6);
  var dateOfMonth = initialDate.slice(6,8);
  var initialDays = ee.Number.parse(ee.Date(startYear.cat('-').cat(month).cat('-').cat(dateOfMonth)).format('DDD'));
  //creating final imgcollection
  var s2_no_cld_shdw =  s2_sr_cld_col_eval//s2_sr_cld_col_eval //s2_sr_col
                        .map(add_cld_shdw_mask)
                        .map(apply_cld_shdw_mask)
                        // .map(apply_scl_water_mask)
                        .map(function(img){
                          return img.clip(boundary)//.divide(10000)
                        })
                        .map(addVariables)
  return s2_no_cld_shdw;                  
}

//methods of process composite image in every 15 days average
var intervalYear= function(start_date,end_date,sampleBounds,preProsS2,intervalValue){
  start_date = ee.Date(start_date)
  end_date = ee.Date(end_date)
  var daysDiff = end_date.difference(start_date, 'days');
  var dayIntervals = ee.List.sequence(0,daysDiff,intervalValue).add(daysDiff).map(function(ele){
                    return ee.Number(ele).toInt();
                    }).distinct();
  var startList = dayIntervals.slice(0,-1);
  var endList = dayIntervals.slice(1,99);
  var pareArrD = ee.Array.cat([startList, endList],1);
  var S2IntervalList = ee.ImageCollection(pareArrD.toList().map(function(pareDay){ 
    var startInterval = start_date.advance(ee.Number(ee.List(pareDay).get(0)),'day');
    var endInterval = start_date.advance(ee.Number(ee.List(pareDay).get(1)),'day');
    var intervalPeriod = ee.Filter.date(startInterval, endInterval);
    var indexSeries = preProsS2.select(['NDVI','GCVI']).filter(intervalPeriod)//'GCVI','NDMI','NDWI','LAI'
    var dayStart=pareArrD.toList().indexOf(ee.List(pareDay)).add(1).multiply(intervalValue)
    var intervalMean = 
    indexSeries.select('NDVI').reduce({reducer:ee.Reducer.mean()}).rename('NDVI')
                            .set('system:day_start',dayStart)
    intervalMean = intervalMean.set('bands',intervalMean.bandNames().length())
    return intervalMean;
  }))
  return S2IntervalList;
}

//incorporating above methods
compositeS2  = function(parameters){
  parameters = parameters || {};
  var startDate = parameters.startDate;
  var endDate = parameters.endDate;
  var S2_CLOUD_FILTER = parameters.S2_CLOUD_FILTER;
  var S2_CLD_PRB_THRESH = parameters.S2_CLD_PRB_THRESH;
  var S2_NIR_DRK_THRESH = parameters.S2_NIR_DRK_THRESH;
  var S2_CLD_PRJ_DIST = parameters.S2_CLD_PRJ_DIST;
  var S2_BUFFER = parameters.S2_BUFFER;
  var area_boundary = parameters.area_boundary;
  var intervalValue = parameters.intervalValue;
  var preProsS2 = S2_cutCldSlw(startDate,endDate,area_boundary,S2_CLOUD_FILTER,S2_CLD_PRB_THRESH,S2_NIR_DRK_THRESH,S2_CLD_PRJ_DIST,S2_BUFFER)
                    .sort('system:day_start');
  var merge15Interval = intervalYear(startDate,endDate,area_boundary,preProsS2,intervalValue);
  return merge15Interval;
}
