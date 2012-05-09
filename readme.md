# Data Science Hackathon (Air Quality Prediction)
This is the code behind our teams entry into the [Kaggle competition](https://www.kaggle.com/c/dsg-hackathon) for this hackathon.

## The team
 * [Ferenc Huszar](mlg.eng.cam.ac.uk/ferenc/)
 * [José Miguel Hernández Lobato](http://www.eng.cam.ac.uk/~jmh233/)
 * [Jedidiah Francis](http://jedifran.com/)
 * [Dan Harvey](https://github.com/danharvey)

## Method
Here there is a short description of what scripts to run and what each script does.

### 1. Estimate missing values and binarize categorical variables.

For this we do

    $ cd dataPreprocessing
    $ R --no-save < preprocessData.R
    $ cd ..

This script loads the training data from "raw_data/TrainingData.csv" and replaces
the categorical variables (day of week, month, hour) by corresponding binary indicators.
The script also estimates the missing values in the data using the variational matrix
factorization method of 

Raiko, T.; Ilin, A. & Juha, K. 
Principal Component Analysis for Large Scale Problems with Lots of Missing Values 
European Conference on Machine Learning: ECML 2007, Springer Berlin / Heidelberg, 2007, 4701, 691-698

The imputed dataset is also standardized so that each non-categorical variable
has zero mean and unit standard deviation.

The input to this script is the file 

 * "raw_data/TrainingData.csv"		->	Contains the original training data.

The output of this script are the files

 * "dataPreprocessing/imputedDataSet.txt"	->	Contains the imputed dataset.
 * "dataPreprocessing/meanMarginals.txt"	-> 	Contains the marginal means of each column from "dataPreprocessing/imputedDataSet.txt".
 * "dataPreprocessing/sdMarginals.txt"	-> 	Contains the marginal standard deviation of each column from "dataPreprocessing/imputedDataSet.txt".

### 2. Reduce the dimensionality of the time-series measurements.

For this we do

    $ cd svd
    $ R --no-save < doSVD.R
    $ cd ..

This script loads the imputed training data from "dataPreprocessing/imputedDataSet.txt" and reduces
the dimensionality of the last 89 columns to just 40 by means of the variational matrix factorization
method described by:

Nakajima, S.; Sugiyama, M. & Tomioka, R.
Global Analytic Solution for Variational Bayesian Matrix Factorization
Advances in Neural Information Processing Systems, 2010

This script generates the file "svd/factors.txt" with the 40 factors that summarize the 89 previous
columns and the file "svd/imputedDataSetWithFactorLoadings.txt" which contains the original imputed
dataset, but with the last 89 columns replaced by the values of the 40 factor loading variables.

The input to this script is the file

 * "dataPreprocessing/imputedDataSet.txt"		->	Contains the imputed dataset.

The output of this script are the files

 * "svd/factors.txt"				->	Contains the 40 factors.
 * "svd/imputedDataSetWithFactorLoadings.txt"	->	Contains the imputed dataset with the 40 factor loading variables.

### 3. Run Python scripts.

    median_hour_month.py

For two chunks there was no test data provided. For these two chunks we couldn't use the method we developed to regress future values against previous observations. These predictions had to be made based on knowledge of the hour and month at which the prediction has to be made.

To handle these two chunks we have computed the global median value of each target variable in each month at certain hours of the day. To make the estimates more accurate we always included data from a month before and after each month. So if we wanted to make a prediction at 10am in October, out output was the median of the target variables at 10 o'clock in September, October and November.

We used median, and not the mean, as we knew the evaluation is going to be based on mean absolute error. The empirical mean would minimise the mean squared error.

The script median_hour_month.py processes the dataset and computes these median values for any month - hour pair.

### 4. Binarize the test file.

For this we do

    $ cd testSetPreprocessing
    $ R --no-save < preprocessData.R
    $ cd ..

This script binarizes the categorical variables in the test file in the same way as we did for the
training file in step 1. For this, the script opens the file "python/training/submission_training.csv" and
generates the file "testSetPreprocessing/binarizedTestDataSet.txt".

The input to this script is the file

"python/training/submission_training.csv"	->	Contains the test set preprocessed by the python scripts so
							that it contains the factor loading variables for each instance.

The output of this script is the file

"testSetPreprocessing/binarizedTestDataSet.txt"	->	The same file as the one above but with the categorical variables binarized.

### 5. Generate a final prediction for the test set.

For this we do

    $ cd fitLinearModel
    $ R --no-save < generateFinalPrediction.R
    $ cd ..

This script fits a linear model that predicts the current values of the factor loading variables given
the previous values of those variables and the current categorical variables
(day of the week, month, hour). These linear models are used to generate the
final predictions for the test set.

The input to this script are the files

 * "imputedDataWithFactorLoadings/factors.txt"		-> 	Contains the 40 factors that sumirze the 89 measurement variables.
 * "testSetPreprocessing/binarizedTestDataSet.txt"		->	Contains the binarized test set.
 * "imputedDataWithFactorLoadings/meanMarginals.txt"	->	Contains the mean of each column of the imputed training set.
 * "imputedDataWithFactorLoadings/sdMarginals.txt"		-> 	Contains the standard deviation of each column of the imputed training set.

 * "python/training/training_*_ubertrain.csv"		-> 	Contains the training data for the prediction task * in advance,
 where * is 1, 2 ,3, 4, 5, 10, 17, 24, 48, and 72.

The output of this script is the file

 * "fitLinearModel/finalResult.txt"			-> 	Contains the final submission file.

### 6. Blending with the benchmark hourly average.

Towards the end of the competition we have reached the point where we thought the predictive power of our simple linear model cannot be improved any further, and we were quite exhausted to try anything more sophisticated. As our predictions were still only marginally better than the provided baseline and sample submission 'submission_hour_chunk_means.csv', we thought we may improve our predictive power by forming an ensemble of our linear model and the provided nonparametric benchmark solution.

blending.py computes a convex combination of our predictions with the supplied sample submission file 'submission_hour_chunk_means.csv' that we downloaded from kaggle.

We have first tried a 95% - 5% blend of our predictions with benchmark predictions. That resulted in a public leaderboard score of 0.22631, thus indeed improving on using our method alone: our solution without bleniding produced a score of 0.22896, and the benchmark solution achieved a score of 0.27532.

We have fitted a parabloa to these three observations (0,0.22896), (0.05, 0.22631) and (1,0.27532), and predicted that if the loss function was quadratic (which we knew it wasn't, but for large number of test predictions it can be reasonably close), an optimal ensemble would blend our method with the benchmark at about 33% - 67% ratio. To avoid overfitting to the part of the test set used for the public leaderboard scores we decided to try a more conservative 25%-75% mixture, which was our final, and best submission with a public score of 0.22075.

The supplied image BlendingInterpolation.pdf shows how our last submission fits well onto a parabola predicted.

