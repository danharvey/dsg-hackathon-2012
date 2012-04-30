
source("vbpcadAnalytic.R")

# We load the data

#dataSet <- as.matrix(read.table("imputedDataSet.txt"))
load("dataSet.dat")

# We take the measurements

measurements <- dataSet[ , 47 : ncol(dataSet) ]

# We do the SVD on these measurements

ret <- fitVBPCADanalytic(measurements, 24)

# We save the factors

write.table(ret$mB, "factors.txt", col.names = F, row.names = F)

# We append the factor loadings to the dataset

dataSet <- cbind(dataSet[ , -seq(47, ncol(dataSet)) ], ret$mA)

# We store the dataset

write.table(dataSet, "imputedDataSetWithFactorLoadings.txt", col.names = F, row.names = F)
