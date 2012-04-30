#
# Author: Jose Miguel Hernandez Lobato
#
# Date: 24 June 2011
#

library(irlba)

##
# Fits the vbpca model given a transaction matrix. The method used
# is the analytical solution described by
#
# S. Nakajima, M. Sugiyama, R. Tomioka, "Global Analytic Solution for Variational Bayesian Matrix Factorization". In NIPS 2010.
#
# @param trans Transaction data as ecoded by the arules package.
# @param H     Dimension of the matrix factorization.
#

fitVBPCADanalytic <- function(trans, H) {

	# We select the dimension of the matrix. The number of rows should be
	# larger than the number of columns.

	if (nrow(trans) > ncol(trans))
		U <- t(trans)
	else
		U <- trans
	L <- nrow(U)
	M <- ncol(U)

	# We compute the svd decomposition of the transaction matrix

	retSVD <- irlba(U, nu = H, nv = H)

	# We initialize sigma to the empirical standard deviation of the SVD residual

	mBaux <- retSVD$u * matrix(retSVD$d, nrow(retSVD$u), ncol(retSVD$u), byrow = T)
	mAaux <- retSVD$v 

	sigma <- sqrt(mean((mBaux %*% t(mAaux) - U)^2))

	sigmaOld <- Inf
	iteration <- 1
#	while (abs(sigma - sigmaOld) > 1e-5) {

		cat(iteration, abs(sigma - sigmaOld), "\n")

		y_h <- retSVD$d
		y_h_under_bar <- (sqrt(L) + sqrt(M)) * sigma

		# We check which components are surely pruned

		nNonZeroComponents <- length(which(y_h >= y_h_under_bar))
		y_h <- y_h[ 1 : nNonZeroComponents ]
		H <- nNonZeroComponents

		# We find the value of c_hat_breve

		c_hat_breve <- sqrt(1 / (2 * L * M) * (y_h^2 - (L + M) * sigma^2 + sqrt((y_h^2 - (L + M) * sigma^2)^2 - 4 * L * M * sigma^4)))

		# We obtain the values of c_a_h and c_b_h

		c_a_h <- sqrt(c_hat_breve)
		c_b_h <- c_a_h

		# We obtain the coefficients of the quartic equation

		e3 <- (L - M)^2 * y_h / (L * M)
		eta_h_squared <- (1 - sigma^2 * L / y_h^2) * (1 - sigma^2 * M / y_h^2) * y_h^2
		e2 <- -(e3 * y_h  + (L^2 + M^2) * eta_h_squared / (L * M) + 2 * sigma^4 / (c_a_h^2 * c_b_h^2))
		e0 <- (eta_h_squared - sigma^4 / (c_a_h^2 * c_b_h^2))^2
		e1 <- e3 * sqrt(e0)

		# We find the value of y_h_hat by solving the quartic equation

		y_h_hat <- rep(0, length(y_h))
		for (i in 1 : length(y_h_hat)) {

			# We solve the quartic equation numerically
		
			roots <- polyroot(c(e0[ i ], e1[ i ], e2[ i ], e3[ i ], 1))

			# We select the second largest real root

			y_h_hat[ i ] <- sort(Re(roots[ abs(Im(roots)) < 1e-6 ]), decreasing = T)[ 2 ]
		}

		# We compute the value of y_h_tilde

		y_h_tilde <- sqrt((L + M) * sigma^2 / 2 + sigma^4 / (2 * c_a_h^2 * c_b_h^2) +
			sqrt(((L + M) * sigma^2 / 2 + sigma^4 / (2 * c_a_h^2 * c_b_h^2))^2 - L * M * sigma^4))

		# We compute the variational Bayes solution

		y_h_hat_VB <- y_h_hat
		y_h_hat_VB[ y_h <= y_h_tilde  ] <- 0

		# We estimate the value of c_h that minimizes the solution

		Delta_h <- M * log(y_h / (M * sigma^2) * y_h_hat_VB + 1) +
			L * log(y_h / (L * sigma^2) * y_h_hat_VB + 1) + 1 / sigma^2 * (-2 * y_h * y_h_hat_VB + L * M * c_hat_breve^2)
		c_h_hat <- c_hat_breve
		c_h_hat[ !(y_h > y_h_under_bar & Delta_h <= 0) ] <- 0

		# We obtain the parameter values of the posterior approximation

		delta_h_hat <- ((M - L) * (y_h - y_h_hat_VB) + sqrt((M - L)^2 * (y_h - y_h_hat_VB)^2 + 4 *
			sigma^4 * L * M / (c_a_h^2 * c_b_h^2))) / (2 * sigma^2 * M * c_a_h^-2)
		eta_h_hat_squared <- eta_h_squared
		eta_h_hat_squared[ y_h <= y_h_tilde ] <- sigma^4 / (c_a_h^2 * c_b_h^2)
		mA <- matrix(sqrt(y_h_hat_VB * delta_h_hat), M, H, byrow = T) * retSVD$v[ , 1 : H ]
		mB <- matrix(sqrt(y_h_hat_VB / delta_h_hat), L, H, byrow = T) * retSVD$u[ , 1 : H ]
		vA <- matrix((-(eta_h_hat_squared - sigma^2 * (M - L)) + sqrt((eta_h_hat_squared - sigma^2 * (M - L))^2 + 4 * M * sigma^2 * eta_h_hat_squared)) /
			(2 * M * (y_h_hat_VB / delta_h_hat + sigma^2 * c_a_h^-2)), M, H, byrow = T)
		vB <- matrix((-(eta_h_hat_squared + sigma^2 * (M - L)) + sqrt((eta_h_hat_squared + sigma^2 * (M - L))^2 + 4 * L * sigma^2 * eta_h_hat_squared)) /
			(2 * L * (y_h_hat_VB * delta_h_hat + sigma^2 * c_b_h^-2)), L, H, byrow = T)

		# We eliminate the columns which are pruned by the ARD method

		H <- length(which(c_h_hat != 0))
		mA <- matrix(mA[ , 1 : H ], nrow(mA), H)
		vA <- matrix(vA[ , 1 : H ], nrow(vA), H)
		mB <- matrix(mB[ , 1 : H ], nrow(mB), H)
		vB <- matrix(vB[ , 1 : H ], nrow(vB), H)
		c_a_h <- c_a_h[ 1 : H ]
		c_b_h <- c_b_h[ 1 : H ]

		# We refine the value of sigma

#		error <- sum(U) + sum(apply(mB^2, 2, sum) * apply(mA^2, 2, sum)) - 2 * sum(mB[index_i, ] * mA[index_j, ])
#		alpha_h <- apply(mA^2, 2, sum) + apply(vA, 2, sum)
#		beta_h <- apply(mB^2, 2, sum) + apply(vB, 2, sum)
#		extraTerm <- sum(alpha_h * beta_h - apply(mA^2, 2, sum) * apply(mB^2, 2, sum))
#		sigmaOld <- sigma
#		sigma <- sqrt(1 / (L * M) * (error + extraTerm))

		iteration <- iteration + 1
#	}

	# We return the model

	if (nrow(trans) > ncol(trans))
		list(mA = mA, vA = vA, mB = mB, vB = vB, sigma = sigma, c_a_h = c_a_h, c_b_h = c_b_h)
	else
		list(mA = mB, vA = vB, mB = mA, vB = vA, sigma = sigma, c_a_h = c_b_h, c_b_h = c_a_h)
}

##
# Computes the score which indicates how likely it is that a

computeVBPCADanalyticScores <- function(set, model) {

	nItems <- nrow(model$mB)
	testTransaction <- rep(0, nItems)
	testTransaction[ set ] <- 1

	H <- ncol(model$mB)

	# We estimate the score

	vPost <- solve(diag(model$c_a_h^-2, H) + model$sigma^-2 * (diag(apply(model$vB, 2, sum), H) + t(model$mB) %*% model$mB))
	mPost <- vPost %*% (t(model$mB) %*% testTransaction) * model$sigma^-2
	as.double(model$mB %*% mPost)
}
