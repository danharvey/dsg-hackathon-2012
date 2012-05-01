
source("vbpcad.R")

# We load the dataset

dataSet <- read.table("../raw_data/TrainingData.csv", header = T, sep = ",")

# We map the categorical variables to binary indicators

dataSetBinaryIndicators <- as.matrix(dataSet[ , 1 : 3 ])

# Indicator function that maps a categorical value to its corresponding binary vector

indicator <- function(x, nValues) { ret <- rep(0, nValues); ret[ x ] <- 1 ; ret }

# We map the month to its category

nValues <- length(unique(dataSet[ , 4 ]))
month <- t(apply(matrix(dataSet[ , 4 ], nrow(dataSet), 1), 1, function(x) indicator(x, nValues)))

dataSetBinaryIndicators <- cbind(dataSetBinaryIndicators, as.matrix(month))

# We map the day of the week to its category

nValues <- length(unique(as.double(dataSet[ , 5 ])))
weekDay <- t(apply(matrix(as.double(dataSet[ , 5 ]), nrow(dataSet), 1), 1, function(x) indicator(x, nValues)))

dataSetBinaryIndicators <- cbind(dataSetBinaryIndicators, as.matrix(weekDay))

# We map the hour to its category

dataSet[ dataSet[ , 6 ] == 0, 6 ] <- 24
nValues <- length(unique(dataSet[ , 6 ]))
hour <- t(apply(matrix(dataSet[ , 6 ], nrow(dataSet), 1), 1, function(x) indicator(x, nValues)))

dataSetBinaryIndicators <- cbind(dataSetBinaryIndicators, as.matrix(hour))

# We keep the column of the column in which the real variables start

startRealVariables <- ncol(dataSetBinaryIndicators) + 1

# We add the rest of the dataset which is not categorical

dataSetBinaryIndicators <- cbind(dataSetBinaryIndicators, as.matrix(dataSet[ , 7 : ncol(dataSet) ]))

# We standardize the marginals to be Gaussian

meanMarginals <- rep(0, ncol(dataSetBinaryIndicators))
sdMarginals <- rep(1, ncol(dataSetBinaryIndicators))
for (i in startRealVariables : ncol(dataSetBinaryIndicators)) {

	# We map the data to the Gaussian domain

	index <- which(!is.na(dataSetBinaryIndicators[ , i ]))
	rawDataMarginal <- dataSetBinaryIndicators[ index, i ]

	meanMarginals[ i ] <- mean(rawDataMarginal)
	sdMarginals[ i ] <- sd(rawDataMarginal)

	rawDataMarginal <- (rawDataMarginal - meanMarginals[ i ]) / sdMarginals[ i ]

	dataSetBinaryIndicators[ index, i ] <- rawDataMarginal
}

# We impute missing data, ignoring the first 3 features: rowID, chunkID, position_within_chunk 

dataSetToImpute <- dataSetBinaryIndicators[ , 4 : ncol(dataSetBinaryIndicators) ]

nNonMissing <- length(dataSetToImpute) - length(which(is.na(dataSetToImpute)))
ratings <- matrix(0, nNonMissing, 3)
counter <- 1
for (i in 1 : nrow(dataSetToImpute)) {
	for (j in 1 : ncol(dataSetToImpute)) {
		if (!is.na(dataSetToImpute[ i, j ])) {
			ratings[ counter, 1 ] <- i
			ratings[ counter, 2 ] <- j
			ratings[ counter, 3 ] <- dataSetToImpute[ i, j ]

			counter <- counter + 1
		}
	}

	print(i)
}

ret <- bvpcadFast(ratings, matrix(0, nrow(dataSetToImpute), 0), matrix(0, ncol(dataSetToImpute), 0), 50)

# We recontstruct the matrix

dataSetImputed <- ret$mPosteriorU %*% t(ret$mPosteriorM)
dataSetImputed <- cbind(dataSetBinaryIndicators[ , 1 : 3 ], dataSetImputed)

index <- which(is.na(dataSetBinaryIndicators))
dataSetBinaryIndicators[ index ] <- dataSetImputed[ index ]

# We store the imputed datset and the means and variances for the marginals

write.table(dataSetBinaryIndicators, "imputedDataSet.txt", col.names = F, row.names = F)

write.table(meanMarginals, "meanMarginals.txt", col.names = F, row.names = F)

write.table(sdMarginals, "sdMarginals.txt", col.names = F, row.names = F)

cat("Done!\n")
