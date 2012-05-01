#########################################################################################################
# Module: epSSreg.R
# Date  : January 2011
# Author: Jose Miguel Hernandez Lobato
# email : josemiguel.hernandez@uam.es
#
# Module that approximates the posterior of the parameters of the the Bayesian
# variable selection model for the linear regression problem. The procedure used to perform
# the approximation is Expectation Propagation.
#
#########################################################################################################
#
# EXPORTED FUNCTIONS: epSSreg
#
#########################################################################################################
#
# The main function that should be used in this module is "epSSreg". You have to call it with the arguments:
#
# 	X -> Design matrix for the regression problem.
#       Y -> Target vector for the regression problem.
#	beta -> Noise precision.
#	p0 -> Prior probability that a feature is relevant for solving the regression problem.
#	v -> Variance of the slab.
#
# "epSSreg" returns the approximate distribution for the posterior as a list with components:
#
#	mw -> Mean vector for the posterior marginal distribution of the model coefficients.
#	vw -> Variances vector for the posterior marginal distribution of the model coefficients.
#	ma -> Mean vector for the posterior marginal distribution of the latent variables.
#	va -> Variances vector for the posterior marginal distribution of the latent variables.
#
#	phi -> Vector with probit transformation of the marginal probabilities of activation of the latent variables.
#	p -> Vector with the marginal prbobabilities of activation of the latent variables.
#
#	t0Hat -> A list with the approximation for the likelihood terms.
#		maHat -> Mean vector of the factorized Gaussian approximation.
#		vaHat -> Variance vector for the factorized Gaussian approximation.
#
#	t1Hat -> A list with the approximation for the delta functions.
#		mwHat -> Mean vector of the factorized Gaussian approximation on the coefficients.
#		vwHat -> Variance vector for the factorized Gaussian approximation on the coefficients.
#		maHat -> Mean vector of the factorized Gaussian approximation on the latent varialbes.
#		vaHat -> Variance vector for the factorized Gaussian approximation on the latent variables.
#
#	t2Hat -> A list with the approximation for the spike and slab prior.
#		mwHat -> Mean vector of the factorized Gaussian approximation on the coefficients.
#		vwHat -> Variance vector for the factorized Gaussian approximation on the coefficients.
#		phiHat -> Probit transformation of the parameter vector for the Bernoulli approximation on the latent variables.
#
#	t3Hat -> A list with the approximation for Bernoulli prior.
#		phiHat -> Probit transformation of the parameter vector for the Bernoulli approximation on the latent variables.
#
#	evidence -> The approximation for the evidence given by EP.
#

v0 <- 0

epSSreg <- function(X, Y, beta = 1, p0 = 0.5, v = 1) {

	d <- ncol(X)
	n <- nrow(X)

	# We initialize the posterior approximation

	t0Hat <- list(maHat = rep(0, n), vaHat = rep(Inf, n))
	t1Hat <- list(mwHat = rep(0, d), vwHat = rep(Inf, d), maHat = rep(0, n), vaHat = rep(Inf, n))
	t2Hat <- list(mwHat = rep(0, d), vwHat = rep(Inf, d), phiHat = rep(0, d))
	t3Hat <- list(phiHat = rep(0, d))
	a <- list(ma = rep(0, d), va = rep(Inf, n), mw = rep(0, d), vw = rep(Inf, d), phi = rep(0, d), p = rep(NA, d),
		t0Hat = t0Hat, t1Hat = t1Hat, t2Hat = t2Hat, t3Hat = t3Hat, indexNegative = c())

	# We process the approximate term for the Bernoulli prior

	a$t3Hat$phiHat <- rep(log(p0 / (1 - p0)), d)
	a$phi <- a$t3Hat$phiHat
	a$p <- rep(p0, d)

	# We process the approximate term for the spike and slab prior

	a$t2Hat$vwHat <- rep(p0 * v, d)
	a$vw <- a$t2Hat$vwHat

	# We process the approximate term for the delta functions

	a$t1Hat$vaHat <- as.double((X^2 * matrix(a$t2Hat$vwHat, n, d, byrow = T)) %*% rep(1, d))
	a$va <- a$t1Hat$vaHat

	# We process the approximate term for the likelihood

	a$t0Hat$maHat <- as.double(Y)
	a$t0Hat$vaHat <- rep(1 / beta, n)
	a$va <- (a$t0Hat$vaHat^-1 + a$t1Hat$vaHat^-1)^-1
	a$ma <- a$va * (a$t0Hat$maHat * a$t0Hat$vaHat^-1 + a$t1Hat$maHat * a$t1Hat$vaHat^-1)

	# Main loop of EP

	i <- 1
	damping <- 0.99
	convergence <- FALSE
	while(!convergence && i < 2) {

		aOld <- a

		# We refine the approximate term for the delta functions

		if (d > n) {
			matrixAux <- (X * matrix(a$t2Hat$vwHat, n, d, byrow = T)) %*% t(X)
			inverseWoodbury <- solve(diag(as.double(a$t0Hat$vaHat)) + matrixAux)
			vectorAux <- a$t2Hat$vwHat^-1 * a$t2Hat$mwHat + t(X) %*% (a$t0Hat$maHat * a$t0Hat$vaHat^-1)
			a$mw <- as.double(a$t2Hat$vwHat * (vectorAux - (t(X) %*% (inverseWoodbury %*% (X %*% (a$t2Hat$vwHat * vectorAux))))))
			a$vw <- as.double(a$t2Hat$vwHat - a$t2Hat$vwHat^2 * (rep(1, n) %*% (X * (inverseWoodbury %*% X))))
			a$ma <- as.double(X %*% a$mw)
			a$va <- as.double(diag(matrixAux - matrixAux %*% inverseWoodbury %*% matrixAux))
		} else {
			Sigma <- chol2inv(chol(diag(as.double(a$t2Hat$vwHat)^-1) + (matrix(a$t0Hat$vaHat^-1, d, n, byrow = T) * t(X)) %*% X))
			a$mw <- as.double(Sigma %*% (a$t2Hat$vwHat^-1 * a$t2Hat$mwHat + t(X) %*% (a$t0Hat$maHat * a$t0Hat$vaHat^-1)))
			a$vw <- as.double(diag(Sigma))
			a$ma <- as.double(X %*% a$mw)
			a$va <- as.double(((X %*% Sigma) * X) %*% rep(1, d))
		}

		a$t1Hat$mwHat <- (damping * (a$mw * a$vw^-1 - a$t2Hat$mwHat * a$t2Hat$vwHat^-1) + (1 - damping) * a$t1Hat$mwHat * a$t1Hat$vwHat^-1)
		a$t1Hat$vwHat <- 1 / (damping * (1 / a$vw - 1 / a$t2Hat$vwHat) + (1 - damping) / a$t1Hat$vwHat)
		a$t1Hat$mwHat <- a$t1Hat$vwHat * a$t1Hat$mwHat

		a$t1Hat$maHat <- (damping * (a$ma * a$va^-1 - a$t0Hat$maHat * a$t0Hat$vaHat^-1) + (1 - damping) * a$t1Hat$maHat * a$t1Hat$vaHat^-1)
		a$t1Hat$vaHat <- 1 / (damping * (1 / a$va - 1 / a$t0Hat$vaHat) + (1 - damping) / a$t1Hat$vaHat)
		a$t1Hat$maHat <- a$t1Hat$vaHat * a$t1Hat$maHat

		# We refine the approximate term for the spike and slab prior

		phiHatNew <- 0.5 * log(a$t1Hat$vwHat + v0) - 0.5 * log(a$t1Hat$vwHat + v) + 0.5 * a$t1Hat$mwHat^2 * ((a$t1Hat$vwHat + v0)^-1 - (a$t1Hat$vwHat + v)^-1)
		aa <- logistic(phiHatNew + a$t3Hat$phiHat) * a$t1Hat$mwHat * (a$t1Hat$vwHat + v)^-1 +
			logistic(-phiHatNew - a$t3Hat$phiHat) * a$t1Hat$mwHat * (a$t1Hat$vwHat + v0)^-1
		bb <- logistic(phiHatNew + a$t3Hat$phiHat) * (a$t1Hat$mwHat^2 - a$t1Hat$vwHat - v) * (a$t1Hat$vwHat + v)^-2 +
			logistic(-phiHatNew - a$t3Hat$phiHat) * (a$t1Hat$mwHat^2 - a$t1Hat$vwHat - v0) * (a$t1Hat$vwHat + v0)^-2

		vwHatNew <- (aa^2 - bb)^-1 - a$t1Hat$vwHat
		mwHatNew <- a$t1Hat$mwHat - aa * (vwHatNew + a$t1Hat$vwHat)

		a$indexNegative <- which(vwHatNew < 0)

			# We minimize the KL divergence with vwHatNew constrained to be positive.

		vwHatNew[ a$indexNegative ] <- 100
		mwHatNew[ a$indexNegative ] <- a$t1Hat$mwHat[ a$indexNegative ] - aa[ a$indexNegative ] *
			(vwHatNew[ a$indexNegative ] + a$t1Hat$vwHat[ a$indexNegative ])

		a$t2Hat$phiHat <- phiHatNew * damping + a$t2Hat$phiHat * (1 - damping)
		a$t2Hat$mwHat <- damping * mwHatNew * vwHatNew^-1 + (1 - damping) * a$t2Hat$mwHat * a$t2Hat$vwHat^-1
		a$t2Hat$vwHat <- 1 / (damping / vwHatNew + (1 - damping) / a$t2Hat$vwHat)
		a$t2Hat$mwHat <- a$t2Hat$mwHat * a$t2Hat$vwHat

			# We update the posterior approximation

		a$vw <- 1 / (1 / a$t1Hat$vwHat + 1 / a$t2Hat$vwHat)
		a$mw <- a$vw * (a$t1Hat$mwHat / a$t1Hat$vwHat + a$t2Hat$mwHat / a$t2Hat$vwHat)
		a$phi <- a$t2Hat$phiHat + a$t3Hat$phiHat
		a$p <- logistic(a$phi)

		# We process the approximate term for the likelihood
		
			# We have to do nothing because for the regression case the approximate term for the likelihood is exact

		# Annealed damping scheme

		damping <- damping * 0.99

		# We check for convergence

		convergence <- checkConvergence(a, aOld)

		i <- i + 1
	}

	# We compute the evidence

#	a$evidence <- computeEvidence(a, Y, X, beta, v)

	# We return the current approximation

	a
}

##
# The logistic function
#

logistic <- function(x) {

	1 / (1 + exp(-x))
}

##
# Checks convergence of the EP algorithm.
#
# Input:
# 	aOld -> The previous approximation.
# 	aNew -> The new approximation.
# Output:
# 	TRUE if the values in aOld are differ from those in aNew by less than a small constant.
#

checkConvergence <- function(aNew, aOld) {

	tol <- 1e-3

	convergence <- max(max(abs(aNew$mw - aOld$mw)))
	convergence <- max(convergence, max(abs(aNew$vw - aOld$vw)))

	print(convergence)

	if (convergence < tol)
		TRUE
	else
		FALSE
}

##
# Function that computes the log evidence
#

computeEvidence <- function(a, Y, X, beta, v) {

	n <- nrow(X)
	d <- ncol(X)

	# We compute the logarithm of s0, s1 and s2

	if (n > d)
		logAlpha <- determinant(t(X) %*% diag(as.double(a$t0Hat$vaHat^-1)) %*% X %*%
			diag(as.double(a$t2Hat$vwHat)) + diag(rep(1, d)))$modulus[[ 1 ]]
	else
		logAlpha <- determinant((matrix(a$t2Hat$vwHat, n, d, byrow = T) * X) %*%
			(t(X) * matrix(a$t0Hat$vaHat^-1, d, n, byrow = T)) + diag(rep(1, n)))$modulus[[ 1 ]]
	
	logs0 <- -n / 2 * log(2 * pi) - 0.5 * sum(log(a$t0Hat$vaHat))

	logs1 <- -n / 2 * log(2 * pi) - 0.5 * sum(log(a$t0Hat$vaHat)) - 0.5 * sum(a$t0Hat$maHat^2 * a$t0Hat$vaHat^-1) -
		0.5 * sum(a$t2Hat$mwHat^2 * a$t2Hat$vwHat^-1) + 0.5 * sum((a$t2Hat$vwHat^-1 * a$t2Hat$mwHat + t(X) %*% (a$t0Hat$maHat * a$t0Hat$vaHat^-1)) * a$mw) -
		0.5 * logAlpha + 1 / 2 * sum(log(1 + a$t2Hat$vwHat * a$t1Hat$vwHat^-1)) + 1 / 2 * sum(log(1 + a$t0Hat$vaHat * a$t1Hat$vaHat^-1)) +
		1 / 2 * sum(a$t2Hat$mwHat^2 * a$t2Hat$vwHat^-1 + a$t1Hat$mwHat^2 * a$t1Hat$vwHat^-1 - a$mw^2 * a$vw^-1) +
		1 / 2 * sum(a$t0Hat$maHat^2 * a$t0Hat$vaHat^-1 + a$t1Hat$maHat^2 * a$t1Hat$vaHat^-1 - a$ma^2 * a$va^-1)

	c <- logistic(a$t3Hat$phiHat) * dnorm(0, a$t1Hat$mwHat, sqrt(a$t1Hat$vwHat + v)) + logistic(-a$t3Hat$phiHat) *
		dnorm(0, a$t1Hat$mwHat, sqrt(a$t1Hat$vwHat))

	logs2 <- sum(log(c) + 1 / 2 * log(1 + a$t1Hat$vwHat * a$t2Hat$vwHat^-1) +
		 1 / 2 * (a$t2Hat$mwHat^2 * a$t2Hat$vwHat^-1 + a$t1Hat$mwHat^2 * a$t1Hat$vwHat^-1 - a$mw^2 * a$vw^-1) +
		 log(logistic(a$phi) / logistic(a$t3Hat$phiHat) + logistic(-a$phi) / logistic(-a$t3Hat$phiHat)))

	aux <- d / 2 * log(2 * pi) + 0.5 * sum(log(a$vw)) - 0.5 * sum(a$t1Hat$mwHat^2 / a$t1Hat$vwHat) -
		0.5 * sum(a$t2Hat$mwHat^2 / a$t2Hat$vwHat) + 0.5 * sum(a$mw^2 / a$vw) +
		n / 2 * log(2 * pi) + 0.5 * sum(log(a$va)) - 0.5 * sum(a$t1Hat$maHat^2 / a$t1Hat$vaHat) -
		0.5 * sum(a$t0Hat$maHat^2 / a$t0Hat$vaHat) + 0.5 * sum(a$ma^2 / a$va)

	logs0 + logs1 + logs2 + aux + sum(log(logistic(a$t2Hat$phiHat) * logistic(a$t3Hat$phiHat) + logistic(-a$t2Hat$phiHat) * logistic(-a$t3Hat$phiHat)))
}
