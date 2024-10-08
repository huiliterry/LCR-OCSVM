/*
* Author: Hui Li
* Institute: CSISS GMU
* Date: August. 2024
* e-mail: hli47@gmu.edu
* This program conducts satellite data process,  OCSVM classifier training, and classification
*/
//import methods as first: compositeS2() and ResCosineCoefs()
//collecting Sentinel-2 images, processing them as composite NDVI series, transforming them to cosine coefficient dataset
var coefficientsImg = function(startDate,endDate,testBounds,S2_CLOUD_FILTER,S2_CLD_PRB_THRESH,S2_NIR_DRK_THRESH,S2_CLD_PRJ_DIST,S2_BUFFER,intervalValue,order){
  var preProsS2 = compositeS2({startDate:startDate,
                                  endDate:endDate,
                                  S2_CLOUD_FILTER:S2_CLOUD_FILTER,
                                  S2_CLD_PRB_THRESH:S2_CLD_PRB_THRESH,
                                  S2_NIR_DRK_THRESH:S2_NIR_DRK_THRESH,
                                  S2_CLD_PRJ_DIST:S2_CLD_PRJ_DIST,
                                  S2_BUFFER:S2_BUFFER,
                                  area_boundary:testBounds,
                                  intervalValue:intervalValue
                  });
  //generating harmonic coefficients
  var coefficients = ResCosineCoefs(order,preProsS2);//Cosine regression coefficients
  return coefficients;
};

//training the OCSVM classifier
var classifierConstruct = function(order,gamma,nu){
  var S2_CLOUD_FILTER = 100;
  var S2_CLD_PRB_THRESH = 20;
  var S2_NIR_DRK_THRESH = 0.25;
  var S2_CLD_PRJ_DIST = 1;
  var S2_BUFFER = 20; 
  var intervalValue =15;
  //preparing the coefficients, training dataset, classifer
  //svm classifier
  var trainClassifier = function(coefficientsImgs,lableCollection){
    // generating the training dataset using the above sample.
    var trainingSample = coefficientsImgs.sampleRegions({
                                            collection: lableCollection,
                                            properties: ['class'],
                                            scale: 10,
                                            tileScale:8
                                          });                           
    // Trianing libsvm classifier
    var classifier = ee.Classifier.libsvm(
        {decisionProcedure:'Margin',svmType:'ONE_CLASS',kernelType:"RBF",gamma:gamma,nu:nu,oneClass:1})
        .train({
          features: trainingSample,
          classProperty: 'class',
          inputProperties: coefficientsImgs.bandNames()
        });
    return classifier;
  };
  // randomly selecting 50 points from each sugarcane field
  var randomSamples = function(feature){
    return ee.FeatureCollection.randomPoints({region: feature.geometry(), points: 50, seed: 0, maxError: 1});
  };
  //burning fields geometry
  var burningSugarcane =  randomSamples(Florida1201_21_1)
                          .merge(randomSamples(Florida1201_21_4))
                          .merge(randomSamples(Florida1126_21_2))
                          .merge(randomSamples(Florida1126_21_5))
                          .merge(randomSamples(Florida1124_21))
                          .merge(randomSamples(Florida1201_21_5))
                          .merge(randomSamples(Florida1126_21_3))
                          .merge(randomSamples(Florida1116_21_1)) 
                          .merge(randomSamples(Florida1116_21_2))
                          .merge(randomSamples(Florida1116_21_3))
                          .merge(randomSamples(Florida1116_21_4))
                          .merge(randomSamples(Florida1116_21_5))
                          .map(function(feature){
                              return feature.set("class",1);
                            }
                          );

  //samples coefficients
  var coefficientsImgs_1 = coefficientsImg('2021-01-19','2021-12-30',Florida1201_21_1.geometry(),S2_CLOUD_FILTER,S2_CLD_PRB_THRESH,S2_NIR_DRK_THRESH,S2_CLD_PRJ_DIST,S2_BUFFER,intervalValue,order)
  var coefficientsImgs_2 = coefficientsImg('2021-01-10','2021-12-30',Florida1201_21_4.geometry(),S2_CLOUD_FILTER,S2_CLD_PRB_THRESH,S2_NIR_DRK_THRESH,S2_CLD_PRJ_DIST,S2_BUFFER,intervalValue,order)
  var coefficientsImgs_3 = coefficientsImg('2021-01-25','2021-12-31',Florida1126_21_2.geometry(),S2_CLOUD_FILTER,S2_CLD_PRB_THRESH,S2_NIR_DRK_THRESH,S2_CLD_PRJ_DIST,S2_BUFFER,intervalValue,order)
  var coefficientsImgs_4 = coefficientsImg('2021-01-05','2021-12-20',Florida1126_21_5.geometry(),S2_CLOUD_FILTER,S2_CLD_PRB_THRESH,S2_NIR_DRK_THRESH,S2_CLD_PRJ_DIST,S2_BUFFER,intervalValue,order)
  var coefficientsImgs_5 = coefficientsImg('2021-01-20','2021-12-31',Florida1124_21.geometry(),S2_CLOUD_FILTER,S2_CLD_PRB_THRESH,S2_NIR_DRK_THRESH,S2_CLD_PRJ_DIST,S2_BUFFER,intervalValue,order)
  var coefficientsImgs_7 = coefficientsImg('2021-01-20','2021-12-31',Florida1201_21_5.geometry(),S2_CLOUD_FILTER,S2_CLD_PRB_THRESH,S2_NIR_DRK_THRESH,S2_CLD_PRJ_DIST,S2_BUFFER,intervalValue,order)
  var coefficientsImgs_8 = coefficientsImg('2021-02-17','2021-12-02',Florida1126_21_3.geometry(),S2_CLOUD_FILTER,S2_CLD_PRB_THRESH,S2_NIR_DRK_THRESH,S2_CLD_PRJ_DIST,S2_BUFFER,intervalValue,order)
  var coefficientsImgs_9 = coefficientsImg('2020-12-05','2021-12-01',Florida1116_21_1.geometry(),S2_CLOUD_FILTER,S2_CLD_PRB_THRESH,S2_NIR_DRK_THRESH,S2_CLD_PRJ_DIST,S2_BUFFER,intervalValue,order)
  var coefficientsImgs_10 = coefficientsImg('2021-01-27','2021-12-12',Florida1116_21_2.geometry(),S2_CLOUD_FILTER,S2_CLD_PRB_THRESH,S2_NIR_DRK_THRESH,S2_CLD_PRJ_DIST,S2_BUFFER,intervalValue,order)
  var coefficientsImgs_12 = coefficientsImg('2021-02-09','2021-12-02',Florida1116_21_4.geometry(),S2_CLOUD_FILTER,S2_CLD_PRB_THRESH,S2_NIR_DRK_THRESH,S2_CLD_PRJ_DIST,S2_BUFFER,intervalValue,order)
  var coefficientsImgs_13 = coefficientsImg('2021-02-04','2021-12-01',Florida1116_21_5.geometry(),S2_CLOUD_FILTER,S2_CLD_PRB_THRESH,S2_NIR_DRK_THRESH,S2_CLD_PRJ_DIST,S2_BUFFER,intervalValue,order)
  var coefficientsImgs = ee.ImageCollection([coefficientsImgs_1,
                                              coefficientsImgs_2,
                                              coefficientsImgs_3,
                                              coefficientsImgs_4,
                                              coefficientsImgs_5,
                                              coefficientsImgs_7,
                                              coefficientsImgs_8,
                                              coefficientsImgs_9,
                                              coefficientsImgs_10,
                                              coefficientsImgs_12,
                                              coefficientsImgs_13,
                                             ]).mosaic();
  /*****************************************/
  print('coefficientsImgs',coefficientsImgs);
  // Map.addLayer(coefficientsImgs,{},'coefficientsImgs')
  //constructing the classifier
  var Classifier = trainClassifier(coefficientsImgs,burningSugarcane);
  // print('Classifier',Classifier)
  return Classifier;
};

//method for sugarcane classification
var classification = function(Classifier,startDate,endDate,areaBounds,S2_CLOUD_FILTER,S2_CLD_PRB_THRESH,S2_NIR_DRK_THRESH,S2_CLD_PRJ_DIST,S2_BUFFER,intervalValue,order){
  testBounds = testBounds.geometry();
  var preProsS2 = compositeS2({startDate:startDate,
                                  endDate:endDate,
                                  S2_CLOUD_FILTER:S2_CLOUD_FILTER,
                                  S2_CLD_PRB_THRESH:S2_CLD_PRB_THRESH,
                                  S2_NIR_DRK_THRESH:S2_NIR_DRK_THRESH,
                                  S2_CLD_PRJ_DIST:S2_CLD_PRJ_DIST,
                                  S2_BUFFER:S2_BUFFER,
                                  area_boundary:testBounds,
                                  intervalValue:intervalValue
                  });
  //generating harmonic coefficients
  var coefficients = ResCoesHarm(order,preProsS2);//Cosine regression
  
  var GCVImax = preProsS2.select('GCVI').reduce(ee.Reducer.max());
  Map.addLayer(GCVImax,{},'GCVImax',false);

  var GCVImaxMask = GCVImax.gte(2.4).and(GCVImax.lte(3.6));
  Map.addLayer(GCVImaxMask,{},'GCVImaxMask',false);
  
  //preparing cropland mask
  var CDL = ee.Image('USDA/NASS/CDL/2022').select("cultivated").clip(testBounds).remap([1,2],[0,1]).selfMask();
  //classification
  var sugarcaneClassification = coefficients.updateMask(CDL)
                                            .classify(Classifier)
                                            .clip(testBounds)
                                            .updateMask(GCVImaxMask)
                                            .unmask()
                                            .toByte();

  //filter and smooth noise
  var kernel = ee.Kernel.square({radius: 1,units:'pixels',normalize:true});
  // Perform an erosion followed by a dilation, display.
  var classificationS = sugarcaneClassification.focalMode({kernel: kernel, iterations: 2}).selfMask().toByte();
  Map.addLayer(classificationS, {min:0,max:1,palette:['black', 'yellow']}, 'classificationS');
  if (classificationS) {
    Export.image.toDrive({
      image: classificationS,
      description: 'fileName',//the description can be re-edit what you want
      folder:'folderName',
      scale: 10,
      crs: 'epsg:5070',
      region: testBounds,
      maxPixels: 10000000000
    });
  }
};
