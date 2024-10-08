/*
* Author: Hui Li
* Institute: CSISS GMU
* Date: August. 2024
* e-mail: hli47@gmu.edu
* This program conducts multi-order cosine regression calculation.
*/
//method of multi-order cosine regression calculation and output coefficient datasets
var ResCosineCoefs = function(order,s2_no_cld_shdw){
  //dependent variable.
  var dependent_NDVI = 'NDVI';  
  // the order of the model
  var cosines = order;
  // Make a list of cosine frequencies to model. 
  var cosineFrequencies = ee.List.sequence(1, cosines);
  // Function to get a sequence of band names for cosine terms.
  var getNames = function(base, list) {
    return ee.List(list).map(function(i) { 
      return ee.String(base).cat(ee.Number(i).int());
    });
  };
  // Construct lists of names for the cosine terms.
  var cosNames = getNames('cos_', cosineFrequencies);
  // Independent variables.
  var independents = ee.List(['constant']).cat(cosNames);
  var addConstant = function(image) {
    return image.addBands(ee.Image(1));
  };
  // Function to add a time band.
  var addTime = function(image) {
    // Compute time in fractional years since the epoch.
    var days = ee.String(image.get('system:day_start'));
    var timeRadians = ee.Image(
      ee.Number(Math.PI).divide(365).multiply(ee.Number.parse(days))
    ).float().rename('t');
    return image.addBands(timeRadians);
  };  
  var addcosines = function(freqs) {
    return function(image) {
      // Make an image of frequencies.
      var frequencies = ee.Image.constant(freqs);
      // This band should represent time in radians.
      var time = ee.Image(image).select('t');
      // var time2 = ee.Image(image).select('tpow2');
      // var ordersT = time.pow(frequencies).rename(tPow);
      // Get the cosine terms.
      var cosines = time.multiply(frequencies).cos().rename(cosNames);
      // Get the sin terms.
      return image.addBands(cosines);
    };
  };
  var cosineS2 = s2_no_cld_shdw
    .map(addConstant)
    .map(addTime)
    .map(addcosines(cosineFrequencies));
  
  var cosineNDVITrend = cosineS2
    .select(independents.add(dependent_NDVI))
    .reduce(ee.Reducer.linearRegression(independents.length(), 1));
    
  // Turn the array image into a multi-band image of coefficients.
  var cosineNDVITrendCoefficients = cosineNDVITrend.select('coefficients')
      .arrayProject([0])
      .arrayFlatten([independents]);
  return cosineNDVITrendCoefficients;
}; 
