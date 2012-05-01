source("epSSreg.R")

# We load the factors

factors <- as.matrix(read.table("../imputedDataWithFactorLoadings/factors.txt"))

# We load the real test data

realTestData <- as.matrix(read.table("../testSetPreprocessing/binarizedTestDataSet.txt"))
realXTest <- realTestData[ , 1 : (43 + 40) ]
#realXTest <- cbind(realXTest, realXTest^2)
YTestPred <- matrix(0, nrow(realXTest), 40)

# We load the marginal means and the marginal standard deviations

meanMarginals <- read.table("../imputedDataWithFactorLoadings/meanMarginals.txt")$V1
sdMarginals <- read.table("../imputedDataWithFactorLoadings/sdMarginals.txt")$V1

error <- c()
counter <- 1
for (j in c(1, 2, 3, 4, 5, 10, 17, 24, 48, 72)) {

	# We load the training data

	trainingData <- as.matrix(read.table(paste("../python/training/training_", j, "_ubertrain.csv", sep = ""), sep = ","))

	X <- trainingData[ , 1 : (43 + 40) ]
#	X <- cbind(X, X^2)
	Y <- trainingData[ , (43 + 40 + 1) : ncol(trainingData) ]

	# We load the test data

#	testData <- as.matrix(read.table(paste("training_", j, "_test.csv", sep = ""), sep = ","))
#	XTest <- testData[ , 1 : (43 + 24) ]
#	XTest <- cbind(XTest, XTest^2)
#	YTest <- testData[ , (ncol(testData) - 39 + 1) : ncol(testData) ]

#	X <- rbind(X, XTest)
#	Y <- rbind(Y, YTest)

	# We fit a fitc model for each target variable

	models <- list()
	for (i in 1 : ncol(Y)) {

		models[[ i ]] <- epSSreg(X, Y[ , i ], beta = 1 / (var(Y[ , i ]) / 4), p0 = 0.9999, v = 1)

		# We make predictions for the real test data

		for (k in seq(counter, nrow(realXTest), 10)) {

			YTestPred[ k , i ] <- realXTest[ k, ] %*% models[[ i ]]$mw
		}

		print(i)
	}

	counter <- counter + 1
}

# We multiply by the factor

YTestDefactored <- YTestPred %*% t(factors)

# We eliminate the normalization

YTestDefactored <- YTestDefactored * matrix(sdMarginals[ 47 : length(sdMarginals) ], nrow(YTestDefactored), ncol(YTestDefactored), byrow = T) +
	matrix(meanMarginals[ 47 : length(meanMarginals) ], nrow(YTestDefactored), ncol(YTestDefactored), byrow = T)

# We select the 39 target variables

YTestDefactored <- YTestDefactored[ , (ncol(YTestDefactored) - 39 + 1) : ncol(YTestDefactored) ]

# We store the results

write.table(YTestDefactored, "finalResult.txt", col.names = F, row.names = F)
