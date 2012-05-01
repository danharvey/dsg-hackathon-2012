/**
 * Module that implements the vbpca diagonal method described by 
 * 
 * Raiko, T.; Ilin, A. & Juha, 
 * K. Kok, J.; Koronacki, J.; Mantaras, R.; Matwin, S.; Mladenic, D. & Skowron, A. (ed.) 
 * Principal Component Analysis for Large Scale Problems with Lots of Missing Values Machine Learning: 
 * ECML 2007, Springer Berlin / Heidelberg, 2007, 4701, 691-698
 *
 * Author: Jose Miguel Hernandez Lobato
 * Date: 2 Apr 2011.
 *
 */

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <stdio.h>
#include <stdlib.h>

#define MAX_PREVIOUS_BOUNDS 10

/* Auxiliary functions */

SEXP getListElement(SEXP list, const char *str);
double evaluateLowerBound(SEXP R_iR, SEXP R_q, SEXP R_dim, double *e, int *eU, int *eM, double *avgTrE, double *avgPV);
void computeGradient(SEXP R_iR, SEXP R_q, SEXP R_dim, double *e, int *eU, int *eM, double *grU, double *grM);
void computeNewPosteriorVariances(SEXP R_iR, SEXP R_q, SEXP R_dim, double *e, int *eU, int *eM, double *vU, double *vM);
void updatePosteriorVariances(SEXP R_ret, SEXP R_dim, double *vU, double *vM);
void updatePosteriorMeans(SEXP R_ret, SEXP R_dim, double *vU, double *vM, double *grU, double *grM, double lrate);
void restorePosteriorMeansAndVariances(SEXP R_ret, SEXP R_dim, double *vU, double *vM, double *grU, double *grM);
void updatePriorAndNoiseVariances(SEXP R_ret, SEXP R_dim, double avgTrE, double avgPV);

double abs2(double x) {

	if (x < 0)
		x = -x;

	return x;
}

/**
 * Main function that implements the optimization process in C 
 *
 * @param iR  List with the information for the ratings.
 * @param q   The initial value of the posterior approximation.
 * @param dim List with the dimensions of the data and the posterior approximation.
 *
 */

SEXP mainOptimizationLoop(SEXP R_iR, SEXP R_q, SEXP R_dim) {

	SEXP R_ret; 
	double bound, boundOld, lrate, avgTrE, avgPV, *p, diffBounds;
	int convergence, iter, i;

	/* We reserve memory for auxiliary data structures */

	double *e = (double *) malloc(sizeof(double) * *INTEGER_POINTER(getListElement(R_dim, "nR")));
	int *eU = (int *) malloc(sizeof(int) * *INTEGER_POINTER(getListElement(R_dim, "nR")));
	int *eM = (int *) malloc(sizeof(int) * *INTEGER_POINTER(getListElement(R_dim, "nR")));

	double *grU = (double *) malloc(sizeof(double) * *INTEGER_POINTER(getListElement(R_dim, "nU")) * *INTEGER_POINTER(getListElement(R_dim, "nF")));
	double *grM = (double *) malloc(sizeof(double) * *INTEGER_POINTER(getListElement(R_dim, "nM")) * *INTEGER_POINTER(getListElement(R_dim, "nF")));

	double *vU = (double *) malloc(sizeof(double) * *INTEGER_POINTER(getListElement(R_dim, "nU")) * *INTEGER_POINTER(getListElement(R_dim, "nF")));
	double *vM = (double *) malloc(sizeof(double) * *INTEGER_POINTER(getListElement(R_dim, "nM")) * *INTEGER_POINTER(getListElement(R_dim, "nF")));

	double *previousBounds = (double *) malloc(sizeof(double) * MAX_PREVIOUS_BOUNDS);

	/* Protection macros for the parameters */

	PROTECT(R_iR = AS_LIST(R_iR));
	PROTECT(R_q = AS_LIST(R_q));
	PROTECT(R_dim = AS_LIST(R_dim));

	/* We create a copy of the current approximation */

	PROTECT(R_ret = duplicate(R_q));

	/* We evaluate the lower bound */

	fprintf(stdout, "Starting the optimization process.\n");
	fflush(stdout);
	boundOld = evaluateLowerBound(R_iR, R_ret, R_dim, e, eU, eM, &avgTrE, &avgPV);

	/* We update the posterior variances */

	computeNewPosteriorVariances(R_iR, R_ret, R_dim, e, eU, eM, vU, vM);
	updatePosteriorVariances(R_ret, R_dim, vU, vM);

	/* We evaluate the lower bound */

	boundOld = evaluateLowerBound(R_iR, R_ret, R_dim, e, eU, eM, &avgTrE, &avgPV);

	/* Main loop of the algorithm */

	lrate = 1e-2;
	iter = 1;
	convergence = 0;
	while (!convergence && iter < 200) {

		/* We compute the new posterior variances and the gradient */

		computeNewPosteriorVariances(R_iR, R_ret, R_dim, e, eU, eM, vU, vM);
		computeGradient(R_iR, R_ret, R_dim, e, eU, eM, grU, grM);

		/* We update the posterior means and then the variances */

		updatePosteriorMeans(R_ret, R_dim, vU, vM, grU, grM, lrate);
		updatePosteriorVariances(R_ret, R_dim, vU, vM);

		/* We evaluate the lower bound */

		bound = evaluateLowerBound(R_iR, R_ret, R_dim, e, eU, eM, &avgTrE, &avgPV);

		if (bound > boundOld) {
			previousBounds[ (iter - 1) %  MAX_PREVIOUS_BOUNDS ] = bound;
			boundOld = bound;
			lrate = lrate * 1.02;

			/* We update the variances for the prior and for the level noise */

			updatePriorAndNoiseVariances(R_ret, R_dim, avgTrE, avgPV); 

		} else {
			previousBounds[ (iter - 1) %  MAX_PREVIOUS_BOUNDS ] = boundOld;
			fprintf(stdout, "Slowing down the learning rate.\n");
			fflush(stdout);
			lrate = lrate * 0.2;
			restorePosteriorMeansAndVariances(R_ret, R_dim, vU, vM, grU, grM);
			bound = evaluateLowerBound(R_iR, R_ret, R_dim, e, eU, eM, &avgTrE, &avgPV);
		}
		
		/* We check for convergence */

		if (iter - 1 >= MAX_PREVIOUS_BOUNDS) {
			diffBounds = 0;
			for (i = 0 ; i < MAX_PREVIOUS_BOUNDS - 1 ; i++) {
				diffBounds += abs2(previousBounds[ i ] - previousBounds[ i + 1 ]);
			}
			diffBounds /= (MAX_PREVIOUS_BOUNDS - 1);
			diffBounds /= abs2(bound);

			if (diffBounds < 1e-7)
				convergence = 1;

			fprintf(stdout, "%4d : %4f : %4f : %6f\n", iter, bound, sqrt(avgTrE), diffBounds); 
			fflush(stdout);
		} else {
			fprintf(stdout, "%4d : %4f : %4f\n", iter, bound, sqrt(avgTrE)); 
			fflush(stdout);
		}

		iter++;
	}

	/* We update the bound value */

	*NUMERIC_POINTER(getListElement(R_ret, "bound")) = bound;

	UNPROTECT(4);

	free(e); free(grU); free(grM); free(eU); free(eM); free(vU); free(vM);

	return R_ret;
}

/* Auxiliary function to get the list element named str, or return NULL */
     
SEXP getListElement(SEXP list, const char *str) {

	SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);

	int i;
     
	for (i = 0; i < length(list); i++)
		if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
			elmt = VECTOR_ELT(list, i);
			break;
	}

	return elmt;
}

/**
 * Function that evaluates the lower bound.
 */

double evaluateLowerBound(SEXP R_iR, SEXP R_q, SEXP R_dim, double *e, int *eU, int *eM, double *avgTrE, double *avgPV) {

	double term1, term2, term3, *mPri, *vPri, *mPosU, *vPosU, *mPosM, *vPosM, *ra, p, v, vNoise;
	int i, j, k, nF, nU, nM, nR, idxU, idxM, *ia, *ja;
	
	nF = *INTEGER_POINTER(getListElement(R_dim, "nF"));
	nU = *INTEGER_POINTER(getListElement(R_dim, "nU"));
	nM = *INTEGER_POINTER(getListElement(R_dim, "nM"));

	/* We compute the value of the first term */

	mPri = NUMERIC_POINTER(getListElement(R_q, "mPriorU"));
	vPri = NUMERIC_POINTER(getListElement(R_q, "vPriorU"));
	mPosU = NUMERIC_POINTER(getListElement(R_q, "mPosteriorU"));
	vPosU = NUMERIC_POINTER(getListElement(R_q, "vPosteriorU"));
	
	term1 = 0;
	for (i = 0 ; i < nU * nF; i++)
		term1 += 0.5 + 0.5 * log(vPosU[ i ] / vPri[ i ]) - (mPosU[ i ] * mPosU[ i ] +
			vPosU[ i ] + mPri[ i ] * mPri[ i ] - 2 * mPri[ i ] * mPosU[ i ]) / (2 * vPri[ i ]);

	/* We compute the value of the second term */

	mPri = NUMERIC_POINTER(getListElement(R_q, "mPriorM"));
	vPri = NUMERIC_POINTER(getListElement(R_q, "vPriorM"));
	mPosM = NUMERIC_POINTER(getListElement(R_q, "mPosteriorM"));
	vPosM = NUMERIC_POINTER(getListElement(R_q, "vPosteriorM"));
	
	term2 = 0;
	for (i = 0 ; i < nM * nF; i++)
		term2 += 0.5 + 0.5 * log(vPosM[ i ] / vPri[ i ]) - (mPosM[ i ] * mPosM[ i ] +
			vPosM[ i ] + mPri[ i ] * mPri[ i ] - 2 * mPri[ i ] * mPosM[ i ]) / (2 * vPri[ i ]);

	/* We compute the value of the third term */

	ia = INTEGER_POINTER(getListElement(R_iR, "csc_ia"));
	ja = INTEGER_POINTER(getListElement(R_iR, "csc_ja"));
	ra = NUMERIC_POINTER(getListElement(R_iR, "csc_ra"));
	vNoise = *NUMERIC_POINTER(getListElement(R_q, "vNoise"));

	/* We run through the movies */

	*avgTrE = 0;
	*avgPV = 0;
	term3 = 0;
	for (i = 0 ; i < nM ; i++) {
		
		/* We run through the users that rated the movie */

		j = ia[ i ] - 1;
		while (j != ia[ i + 1 ] - 1) {

			/* We compute the prediction and the variance of the prediction and update the value of the third term */

			p = v = 0;
			for (k = 0 ; k < nF ; k++) {
				idxU = ja[ j ] - 1 + k * nU;
				idxM = i + k * nM;
				p += mPosU[ idxU ] * mPosM[ idxM ];
				v += mPosU[ idxU ] * mPosU[ idxU ] * vPosM[ idxM ] +
					mPosM[ idxM ] * mPosM[ idxM ] * vPosU[ idxU ] +
					vPosM[ idxM ] * vPosU[ idxU ];
			}

			/* We store the error value and update the value of the third term */
	
			e[ j ] = ra[ j ] - p; eU[ j ] = ja[ j ] - 1; eM[ j ] = i;
			term3 += (e[ j ] * e[ j ] + v) / (2 * vNoise) + 0.5 * log(2 * PI * vNoise);

			/* We update the average training error and the average predictive variance */
		
			*avgTrE += e[ j ] * e[ j ];
			*avgPV += v;

			j++;

		} /* while */

	} /* for */

	/* We normalize the sum of errors and predictive variances */

	nR = *INTEGER_POINTER(getListElement(R_dim, "nR"));

	*avgTrE /= nR;
	*avgPV /= nR;

	return term1 + term2 - term3;
}

/**
 * Function that computes the gradient of lower bound with respect to the posterior means.
 */

void computeGradient(SEXP R_iR, SEXP R_q, SEXP R_dim, double *e, int *eU, int *eM, double *grU, double *grM) {

	double *mPri, *vPri, *mPosU, *vPosU, *mPosM, *vPosM, vNoise;
	int i, j, nF, nU, nM, nR;

	nF = *INTEGER_POINTER(getListElement(R_dim, "nF"));
	nU = *INTEGER_POINTER(getListElement(R_dim, "nU"));
	nM = *INTEGER_POINTER(getListElement(R_dim, "nM"));
	nR = *INTEGER_POINTER(getListElement(R_dim, "nR"));

	mPosU = NUMERIC_POINTER(getListElement(R_q, "mPosteriorU"));
	vPosU = NUMERIC_POINTER(getListElement(R_q, "vPosteriorU"));
	mPosM = NUMERIC_POINTER(getListElement(R_q, "mPosteriorM"));
	vPosM = NUMERIC_POINTER(getListElement(R_q, "vPosteriorM"));

	/* We compute the gradient for the users */

	mPri = NUMERIC_POINTER(getListElement(R_q, "mPriorU"));
	vPri = NUMERIC_POINTER(getListElement(R_q, "vPriorU"));

	for (i = 0 ; i < nU * nF; i++)
		grU[ i ] = (mPri[ i ] - mPosU[ i ]) / vPri[ i ];

	vNoise = *NUMERIC_POINTER(getListElement(R_q, "vNoise"));

	/* We run through the errors and the features */

	for (i = 0 ; i < nR ; i++) {
		for (j = 0 ; j < nF ; j++) {
			grU[ eU[ i ] + j * nU ] += (e[ i ] * mPosM[ eM[ i ] + j * nM ] - vPosM[ eM[ i ] + j * nM ] * mPosU[ eU[ i ] + j * nU ]) / vNoise;
		}
	}

	/* We compute the gradient for the movies */

	mPri = NUMERIC_POINTER(getListElement(R_q, "mPriorM"));
	vPri = NUMERIC_POINTER(getListElement(R_q, "vPriorM"));

	for (i = 0 ; i < nM * nF; i++)
		grM[ i ] = (mPri[ i ] - mPosM[ i ]) / vPri[ i ];

	vNoise = *NUMERIC_POINTER(getListElement(R_q, "vNoise"));

	/* We run through the errors and the features */

	for (i = 0 ; i < nR ; i++) {
		for (j = 0 ; j < nF ; j++) {
			grM[ eM[ i ] + j * nM ] += (e[ i ] * mPosU[ eU[ i ] + j * nU ] - vPosU[ eU[ i ] + j * nU ] * mPosM[ eM[ i ] + j * nM ]) / vNoise;
		}
	}
}

/**
 * Function that computes the update for the posterior variances.
 */

void computeNewPosteriorVariances(SEXP R_iR, SEXP R_q, SEXP R_dim, double *e, int *eU, int *eM, double *vU, double *vM) {

	double *mPri, *vPri, *mPosU, *vPosU, *mPosM, *vPosM, vNoise;
	int i, j, nF, nU, nM, nR;

	nF = *INTEGER_POINTER(getListElement(R_dim, "nF"));
	nU = *INTEGER_POINTER(getListElement(R_dim, "nU"));
	nM = *INTEGER_POINTER(getListElement(R_dim, "nM"));
	nR = *INTEGER_POINTER(getListElement(R_dim, "nR"));

	mPosU = NUMERIC_POINTER(getListElement(R_q, "mPosteriorU"));
	vPosU = NUMERIC_POINTER(getListElement(R_q, "vPosteriorU"));
	mPosM = NUMERIC_POINTER(getListElement(R_q, "mPosteriorM"));
	vPosM = NUMERIC_POINTER(getListElement(R_q, "vPosteriorM"));

	/* We compute the new variances for the users */

	vPri = NUMERIC_POINTER(getListElement(R_q, "vPriorU"));

	for (i = 0 ; i < nU * nF; i++)
		vU[ i ] = 1 / vPri[ i ];

	vNoise = *NUMERIC_POINTER(getListElement(R_q, "vNoise"));

	/* We run through the errors and the features */

	for (i = 0 ; i < nR ; i++)
		for (j = 0 ; j < nF ; j++)
			vU[ eU[ i ] + j * nU ] += (mPosM[ eM[ i ] + j * nM ] * mPosM[ eM[ i ] + j * nM ] + vPosM[ eM[ i ] + j * nM ]) / vNoise;

	/* We compute the inverse of the variance */

	for (i = 0 ; i < nU * nF; i++)
		vU[ i ] = 1 / vU[ i ];

	/* We compute the new variances for the movies */

	vPri = NUMERIC_POINTER(getListElement(R_q, "vPriorM"));

	for (i = 0 ; i < nM * nF; i++)
		vM[ i ] = 1 / vPri[ i ];

	/* We run through the errors and the features */

	for (i = 0 ; i < nR ; i++)
		for (j = 0 ; j < nF ; j++)
			vM[ eM[ i ] + j * nM ] += (mPosU[ eU[ i ] + j * nU ] * mPosU[ eU[ i ] + j * nU ] + vPosU[ eU[ i ] + j * nU ]) / vNoise;

	/* We compute the inverse of the variance */

	for (i = 0 ; i < nM * nF; i++)
		vM[ i ] = 1 / vM[ i ];

	/* The process is repeated because one update depends on the other */

	vPri = NUMERIC_POINTER(getListElement(R_q, "vPriorU"));
	for (i = 0 ; i < nU * nF; i++)
		vU[ i ] = 1 / vPri[ i ];
	for (i = 0 ; i < nR ; i++)
		for (j = 0 ; j < nF ; j++)
			vU[ eU[ i ] + j * nU ] += (mPosM[ eM[ i ] + j * nM ] * mPosM[ eM[ i ] + j * nM ] + vM[ eM[ i ] + j * nM ]) / vNoise;
	for (i = 0 ; i < nU * nF; i++)
		vU[ i ] = 1 / vU[ i ];

	vPri = NUMERIC_POINTER(getListElement(R_q, "vPriorM"));
	for (i = 0 ; i < nM * nF; i++)
		vM[ i ] = 1 / vPri[ i ];
	for (i = 0 ; i < nR ; i++)
		for (j = 0 ; j < nF ; j++)
			vM[ eM[ i ] + j * nM ] += (mPosU[ eU[ i ] + j * nU ] * mPosU[ eU[ i ] + j * nU ] + vU[ eU[ i ] + j * nU ]) / vNoise;
	for (i = 0 ; i < nM * nF; i++)
		vM[ i ] = 1 / vM[ i ];
}

/**
 * Function that updates the value of the posterior variances.
 */

void updatePosteriorVariances(SEXP R_ret, SEXP R_dim, double *vU, double *vM) {

	double *vPosU, *vPosM, backup;
	int i, nF, nU, nM;

	nF = *INTEGER_POINTER(getListElement(R_dim, "nF"));
	nU = *INTEGER_POINTER(getListElement(R_dim, "nU"));
	nM = *INTEGER_POINTER(getListElement(R_dim, "nM"));

	vPosU = NUMERIC_POINTER(getListElement(R_ret, "vPosteriorU"));
	vPosM = NUMERIC_POINTER(getListElement(R_ret, "vPosteriorM"));

	/* We update the posterior variances */

	for (i = 0 ; i < nU * nF; i++) {
		backup = vPosU[ i ];
		vPosU[ i ] = vU[ i ];
		vU[ i ] = backup;
	}
	for (i = 0 ; i < nM * nF; i++) {
		backup = vPosM[ i ];
		vPosM[ i ] = vM[ i ];
		vM[ i ] = backup;
	}
}

/**
 * Function that updates the posterior means.
 */

void updatePosteriorMeans(SEXP R_ret, SEXP R_dim, double *vU, double *vM, double *grU, double *grM, double lrate) {

	double *mPosU, *mPosM, backup, alpha = 2.0 / 3;
	int i, nF, nU, nM;

	nF = *INTEGER_POINTER(getListElement(R_dim, "nF"));
	nU = *INTEGER_POINTER(getListElement(R_dim, "nU"));
	nM = *INTEGER_POINTER(getListElement(R_dim, "nM"));

	mPosU = NUMERIC_POINTER(getListElement(R_ret, "mPosteriorU"));
	mPosM = NUMERIC_POINTER(getListElement(R_ret, "mPosteriorM"));

	/* We update the posterior means */

	for (i = 0 ; i < nU * nF ; i++) {
		backup = mPosU[ i ];
		mPosU[ i ] += lrate * pow(vU[ i ], alpha) * grU[ i ];
		grU[ i ] = backup;
	}

	for (i = 0 ; i < nM * nF ; i++) {
		backup = mPosM[ i ];
		mPosM[ i ] += lrate * pow(vM[ i ], alpha) * grM[ i ];
		grM[ i ] = backup;
	}
}

/**
 * Function that restores the psterior means and variances.
 */

void restorePosteriorMeansAndVariances(SEXP R_ret, SEXP R_dim, double *vU, double *vM, double *grU, double *grM) {

	double *mPosU, *mPosM, *vPosU, *vPosM;
	int i, nF, nU, nM;

	nF = *INTEGER_POINTER(getListElement(R_dim, "nF"));
	nU = *INTEGER_POINTER(getListElement(R_dim, "nU"));
	nM = *INTEGER_POINTER(getListElement(R_dim, "nM"));

	mPosU = NUMERIC_POINTER(getListElement(R_ret, "mPosteriorU"));
	mPosM = NUMERIC_POINTER(getListElement(R_ret, "mPosteriorM"));
	vPosU = NUMERIC_POINTER(getListElement(R_ret, "vPosteriorU"));
	vPosM = NUMERIC_POINTER(getListElement(R_ret, "vPosteriorM"));

	/* We restore the posterior means and variances */

	for (i = 0 ; i < nU * nF ; i++) {
		mPosU[ i ] = grU[ i ];
		vPosU[ i ] = vU[ i ];
	}

	for (i = 0 ; i < nM * nF ; i++) {
		mPosM[ i ] = grM[ i ];
		vPosM[ i ] = vM[ i ];
	}
}

/**
 * Function updates the variances of the prior and the level of noise.
 */

void updatePriorAndNoiseVariances(SEXP R_ret, SEXP R_dim, double avgTrE, double avgPV) {

	double *vNoise, *mPosU, *mPosM, *vPosU, *vPosM, *mPriU, *mPriM, *vPriU,*vPriM, newVariance;
	int i, j, nU, nM, nUF, nMF, nLF, nF;

	/* We update the variance for the level of noise and the variances of the prior if necessary */

	if (*INTEGER_POINTER(getListElement(R_ret, "optimizeNoiseAndPrior")) == 1) {

		vNoise = NUMERIC_POINTER(getListElement(R_ret, "vNoise"));
		*vNoise = avgTrE + avgPV;

		mPosU = NUMERIC_POINTER(getListElement(R_ret, "mPosteriorU"));
		mPosM = NUMERIC_POINTER(getListElement(R_ret, "mPosteriorM"));
		vPosU = NUMERIC_POINTER(getListElement(R_ret, "vPosteriorU"));
		vPosM = NUMERIC_POINTER(getListElement(R_ret, "vPosteriorM"));
		mPriU = NUMERIC_POINTER(getListElement(R_ret, "mPriorU"));
		mPriM = NUMERIC_POINTER(getListElement(R_ret, "mPriorM"));
		vPriU = NUMERIC_POINTER(getListElement(R_ret, "vPriorU"));
		vPriM = NUMERIC_POINTER(getListElement(R_ret, "vPriorM"));

		nU  = *INTEGER_POINTER(getListElement(R_dim, "nU"));
		nM  = *INTEGER_POINTER(getListElement(R_dim, "nM"));
		nF  = *INTEGER_POINTER(getListElement(R_dim, "nF"));
		nUF = *INTEGER_POINTER(getListElement(R_dim, "nUF"));
		nMF = *INTEGER_POINTER(getListElement(R_dim, "nMF"));
		nLF = *INTEGER_POINTER(getListElement(R_dim, "nLF"));
	
		/* We update the variance of the prior for the coefficients of the users that match the real features */
		
		for (i = 0 ; i < nMF ; i++) {
	
			newVariance = 0;
			for (j = 0 ; j < nU ; j++)
				newVariance += vPosU[ j + i * nU ] + (mPosU[ j + i * nU ] - mPriU[ j + i * nU ]) * (mPosU[ j + i * nU ] - mPriU[ j + i * nU ]);
			newVariance /= nU;
	
			for (j = 0 ; j < nU ; j++)
				vPriU[ j + i * nU ] = newVariance;
		}
	
		/* We update the variance of the prior for the latent features of the users */
	
		for (i = nMF + nUF ; i < nF ; i++) {
	
			newVariance = 0;
			for (j = 0 ; j < nU ; j++)
				newVariance += vPosU[ j + i * nU ] + (mPosU[ j + i * nU ] - mPriU[ j + i * nU ]) * (mPosU[ j + i * nU ] - mPriU[ j + i * nU ]);
			newVariance /= nU;
	
			for (j = 0 ; j < nU ; j++)
				vPriU[ j + i * nU ] = newVariance;
		}
	
		/* We update the variance of the prior for the coefficients of the movies that match the real features */
		
		for (i = nMF ; i < nMF + nUF ; i++) {
	
			newVariance = 0;
			for (j = 0 ; j < nM ; j++)
				newVariance += vPosM[ j + i * nM ] + (mPosM[ j + i * nM ] - mPriM[ j + i * nM ]) * (mPosM[ j + i * nM ] - mPriM[ j + i * nM ]);
			newVariance /= nM;
	
			for (j = 0 ; j < nM ; j++)
				vPriM[ j + i * nM ] = newVariance;
		}
	}
}
