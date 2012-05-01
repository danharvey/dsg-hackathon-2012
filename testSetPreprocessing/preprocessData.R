
# We load the dataset

dataSet <- as.matrix(read.table("../python/training/submission_training.csv", sep = ","))

# We map the categorical variables to binary indicators

dataSetBinaryIndicators <- matrix(0, nrow(dataSet), 0)

# Indicator function that maps a categorical value to its corresponding binary vector

indicator <- function(x, nValues) { ret <- rep(0, nValues); ret[ x ] <- 1 ; ret }

# We map the month to its category

nValues <- 12
month <- t(apply(matrix(as.double(dataSet[ , 1 ]), nrow(dataSet), 1), 1, function(x) indicator(x, nValues)))

dataSetBinaryIndicators <- cbind(dataSetBinaryIndicators, as.matrix(month))

# We map the hour to its category

dataSet[ dataSet[ , 4 ] == 0, 6 ] <- 24
nValues <- 24
hour <- t(apply(matrix(as.double(dataSet[ , 2 ]), nrow(dataSet), 1), 1, function(x) indicator(x, nValues)))

dataSetBinaryIndicators <- cbind(dataSetBinaryIndicators, as.matrix(hour))

# We map the day of the week to its category

nValues <- 7
weekDay <- t(apply(matrix(as.double(as.factor(dataSet[ , 3 ])), nrow(dataSet), 1), 1, function(x) indicator(x, nValues)))

dataSetBinaryIndicators <- cbind(dataSetBinaryIndicators, as.matrix(weekDay))

# We add the rest of the dataset which is not categorical

dataSetBinaryIndicators <- cbind(dataSetBinaryIndicators, as.matrix(dataSet[ , 4 : ncol(dataSet) ]))

dataSetBinaryIndicators <- matrix(as.double(dataSetBinaryIndicators), nrow(dataSetBinaryIndicators), ncol(dataSetBinaryIndicators))

write.table(dataSetBinaryIndicators, "binarizedTestDataSet.txt", col.names = F, row.names = F)

cat("Done!\n")
