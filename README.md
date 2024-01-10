# LCR-OCSVM
Sugarcane transfer learning classification using Linear Cosine Regression and One-Class Support Vector Machine

The workflow starts at Sentinel-2 Image Processing that collects image during a specific period, removes cloud and corresponding shadow, and generates composite images.
The core program is "Calculating Linear Regression", fitting NDVI series, extracting regression coefficients.
The burning sugarcane fields were utilized by training sample process in "Training and Classification". The OCSVM is contrusted, conducting the classification through spatiotemporal.
