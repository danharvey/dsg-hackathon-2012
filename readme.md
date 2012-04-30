# Data Science Hackathon (Air Quality Prediction)

This is the code behind our teams entry into the [Kaggle competition](https://www.kaggle.com/c/dsg-hackathon) for this hackathon.

## The team
 * [Ferenc Huszar](mlg.eng.cam.ac.uk/ferenc/)
 * [José Miguel Hernández Lobato](http://www.eng.cam.ac.uk/~jmh233/)
 * [Jedidiah Francis](http://jedifran.com/)
 * [Dan Harvey](https://github.com/danharvey)

## Method
 - Imputed data and reduce factors with SVD - [Principal Component Analysis for Large Scale Problems with Lots of Missing Values](http://www.mendeley.com/research/principal-component-analysis-for-large-scale-problems-with-lots-of-missing-values/)

 - Linear (least squares?) predictors to SVD factors
   - We trained one predictor per hour range ahead per factor.
   - Base on previous hours data only

 - Created train/test set for each of these range ahead
   - Created as many example as possible inside the hours of the 8 days given
   - Split the train:test at 2:1 randomly assigning the examples.

 - Baseline from the example code for month / hour averages

 - Blending in the time-series prediction and baseline using a quadratic combination
 
 - Optimise the weights on the quadratic combination using our test set

## Create submission

To create a submission in the for the form of sending to Kaggle you do the following 

    make
    
To get clean up everything this produced run

    make clean