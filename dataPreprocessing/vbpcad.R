
library(Matrix)

##
# Module that implements the vbpca diagonal method described by 
#
# Raiko, T.; Ilin, A. & Juha, 
# K. Kok, J.; Koronacki, J.; Mantaras, R.; Matwin, S.; Mladenic, D. & Skowron, A. (ed.) 
# Principal Component Analysis for Large Scale Problems with Lots of Missing Values Machine Learning: 
# ECML 2007, Springer Berlin / Heidelberg, 2007, 4701, 691-698
#
# Author: Jose Miguel Hernandez Lobato
# Date: 2 Apr 2011.
#
# The main function exported by this module is
#
# 	bvpcadFast <- function(r, uF, mF, nLF)
#
# which admits:
#
#	r   -> A matrix of ratings with rows (user_id, movie_id, rating_value)0e
#	uF  -> Features for the users, one row per user, one column per feature.
#	mF  -> Features for the movies, one row per movie, one column per feature.
#	nLF -> Numer of latent features to introduce into the model.
#
# The function returns a posterior approximation with the following fields:
#
#	mPosteriorU -> matrix with the posterior mean values for the users. One row per user, one column per feature.
#	mPosteriorM -> matrix with the posterior mean values for the movies. One row per movie, one column per feature.
#	vPosteriorU -> matrix with the posterior variance values for the users. One row per user, one column per feature.
#	vPosteriorM -> matrix with the posterior variance values for the movies. One row per movie, one column per feature.
#
#	mPriorU -> matrix with the prior mean values for the users. One row per user, one column per feature.
#	mPriorM -> matrix with the prior mean values for the movies. One row per movie, one column per feature.
#	vPriorU -> matrix with the prior variance values for the users. One row per user, one column per feature.
#	vPriorM -> matrix with the prior variance values for the movies. One row per movie, one column per feature.
#
 
library(SparseM)

library(irlba)

dyn.load("vbpcad.so")

##
# Function that initializes the current posterior approximation
#

initializeQ <- function(realUserFeatures, realMovieFeatures, nLatentFeatures, ratingsM, optimizeNoiseAndPrior, priorInitialization) {

	cat("Initializing the posterior approximation randomly.\n")

	nUsers <- nrow(realUserFeatures)
	nMovies <- nrow(realMovieFeatures)
	nRealUserFeatures <- ncol(realUserFeatures)
	nRealMovieFeatures <- ncol(realMovieFeatures)

	# We initialize the prior means and variances
	#
	#	mPriorUsers:  (coefficients for observed movie features, observed user features, extra latent features)
	#	mPriorMovies: (observed movie features, coefficients for observed user features, extra latent features)
	#

	if (is.null(priorInitialization)) {
		mPriorU <- cbind(matrix(0, nUsers, nRealMovieFeatures), realUserFeatures, matrix(0, nUsers, nLatentFeatures))
		mPriorM <- cbind(realMovieFeatures, matrix(0, nMovies, nRealUserFeatures), matrix(0, nMovies, nLatentFeatures))

		vPriorU <- cbind(matrix(1, nUsers, nRealMovieFeatures), matrix(1e-6, nUsers, nRealUserFeatures), matrix(1, nUsers, nLatentFeatures))
		vPriorM <- cbind(matrix(1e-6, nMovies, nRealMovieFeatures), matrix(1, nMovies, nRealUserFeatures),  matrix(1, nMovies, nLatentFeatures))
	} else {
		mPriorU <- priorInitialization$mPriorU
		vPriorU <- priorInitialization$vPriorU
		mPriorM <- priorInitialization$mPriorM
		vPriorM <- priorInitialization$vPriorM
	}

	# We initialize the posterior means
	# If the number of real features is 0 we use the svd method

	if (nRealUserFeatures == 0 && nRealUserFeatures == 0) {

		ratingsM <- as.matrix.coo(ratingsM)
		ratingsM <- sparseMatrix(i = ratingsM@ia, j = ratingsM@ja, x = ratingsM@ra, dims = ratingsM@dimension)
		ret <- irlba(ratingsM, nu = nLatentFeatures, nv = nLatentFeatures)
		mPosteriorU <- ret$u %*% diag(ret$d, length(ret$d))
		mPosteriorM <- ret$v

	} else {

		mPosteriorU <- cbind(matrix(rnorm(nUsers * nRealMovieFeatures), nUsers, nRealMovieFeatures), realUserFeatures,
			matrix(rnorm(nUsers * nLatentFeatures), nUsers, nLatentFeatures))

		mPosteriorM <- cbind(realMovieFeatures, matrix(rnorm(nMovies * nRealUserFeatures), nMovies, nRealUserFeatures),
		matrix(rnorm(nMovies * nLatentFeatures), nMovies, nLatentFeatures))
	}

	# We initialize the posterior variances

	vPosteriorU <- cbind(matrix(1, nUsers, nRealMovieFeatures), matrix(1e-6, nUsers, nRealUserFeatures), matrix(1, nUsers, nLatentFeatures))
	vPosteriorM <- cbind(matrix(1e-6, nMovies, nRealMovieFeatures), matrix(1, nMovies, nRealUserFeatures),  matrix(1, nMovies, nLatentFeatures))

	# We initialize the variance of the additive noise

	vNoise <- 1

	# We return the posterior approximation

	list(mPriorU = mPriorU, vPriorU = vPriorU, mPriorM = mPriorM, vPriorM = vPriorM, mPosteriorU = mPosteriorU, mPosteriorM = mPosteriorM,
		vPosteriorU = vPosteriorU, vPosteriorM = vPosteriorM, vNoise = vNoise, bound = -Inf, optimizeNoiseAndPrior = as.integer(optimizeNoiseAndPrior))
}

##
# The VBPCAd method, calling the C routine for the optimization process.
#
# @param r   Matrix with the ratings, one row per rating. Rating format (user_id, movie_id, rating).
# @param uF  Matrix with the real features for the users.
# @param uM  Matrix with the real features for the movies.
# @param nLF Number of latent features to use.
# @param q   Current posterior approximation. NULL if there is no one.
#

bvpcadFast <- function(r, uF, mF, nLF, q = NULL, optimizeNoiseAndPrior = 1, priorInitialization = NULL) {

	# We extend the real features with a bias

#	uF <- cbind(uF, rep(1, nrow(uF)))
#	mF <- cbind(mF, rep(1, nrow(mF)))

	nUF  <- ncol(uF)
	nMF  <- ncol(mF)

	# We store the ratings in a sparse matrix with dimension nUsers x nMovies
	
	rM <- as.matrix.csr(new("matrix.coo", ia = as.integer(r[ , 1 ]), ja = as.integer(r[ , 2 ]), ra = r[ , 3 ], dimension = c(nrow(uF), nrow(mF))))

	# We initialize the posterior approximation

	if (is.null(q))
		q <- initializeQ(uF, mF, nLF, rM, optimizeNoiseAndPrior, priorInitialization)

	# We do the optimization

	optimizeFast(rM, q, nUF, nMF, nLF)
}

##
# Function that calls the C implementation of the main loop of the optimization process.
#
# @param rM  Matrix with the ratings in sparse csr format. As many rows as users.
# @param q   Current posterior approximation.
# @param nUF Number of real features for the users.
# @param uMF Number of real features for the movies.
# @param nLF Number of latent features.
#

optimizeFast <- function(rM, q, nUF, nMF, nLF) {

	# We initialize the parameters to pass to the C optimizer

	dim <- list(nLF = as.integer(nLF), nUF = nUF, nMF = nMF, nM = ncol(rM), nU = nrow(rM), nF = as.integer(nUF + nMF + nLF), nR = length(rM@ra))

	rMcsc <- as.matrix.csc(rM)
	iR <- list(csr_ra = rM@ra, csr_ja = rM@ja, csr_ia = rM@ia, csc_ra = rMcsc@ra, csc_ja = rMcsc@ja, csc_ia = rMcsc@ia)

	ret <- .Call("mainOptimizationLoop", iR, q, dim)

	ret
}
