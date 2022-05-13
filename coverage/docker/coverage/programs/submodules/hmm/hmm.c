#include <math.h>
#include <assert.h>
#include "data_types.h"
#include "hmm.h"
#include <stdio.h>
#include "block_it.h"
#include "common.h"
#include "sonLib.h"

#define PI 3.14159

static pthread_mutex_t chunkMutex;

/**
 * Constructs a Gaussian object containing the mean vectors, covariance
 * matrices and the weights for all of its mixture components
 *
 * @param mu	An array of double vectors containing the initial mean vectors
 * @param cov	An array of double matrices containing the initial covariance matrices
 * @param n	The number of mixture components 
 * @return gaussian
 */
Gaussian* Gaussian_construct(VectorDouble** mu, MatrixDouble** cov, int n){
	Gaussian* gaussian = malloc(1 * sizeof(Gaussian));
	// Mixture Gaussian Parameters:
	//
	// Allocating arrays for gaussian parameters
	gaussian->mu = VectorDouble_constructArray1D(n, mu[0]->dim);
	gaussian->cov = MatrixDouble_constructArray1D(n, mu[0]->dim, mu[0]->dim);
	// Allocating and initializing mixture weights
	gaussian->weights = malloc(n * sizeof(double));
	gaussian->n = n;
	for (int m = 0; m < n; m++){
		// init mean vector
		VectorDouble_copyInPlace(gaussian->mu[m], mu[m]);
		// init covariance matrix
		MatrixDouble_copyInPlace(gaussian->cov[m], cov[m]);
		// init weights uniformly
		gaussian->weights[m] = 1.0 / n;
	}
	// Estimation Statistics:
	//
	// Allocating and initializing arrays for estimating parameters
        gaussian->muNum = VectorDouble_constructArray1D(n, mu[0]->dim);
        gaussian->muDenom = VectorDouble_constructArray1D(n, mu[0]->dim);
        gaussian->covNum = MatrixDouble_constructArray1D(n, mu[0]->dim, mu[0]->dim);
        gaussian->covDenom = (double*) malloc(n * sizeof(double));
        memset(gaussian->covDenom, 0, n * sizeof(double));
	gaussian->weightNum = (double*) malloc(n * sizeof(double));
        memset(gaussian->weightNum, 0, n * sizeof(double));
	return gaussian;
}

/**
 * Destructs a Gaussian object
 *
 * @param gaussian
 */
void Gaussian_destruct(Gaussian* gaussian){
        VectorDouble_destructArray1D(gaussian->mu, gaussian->n);
        VectorDouble_destructArray1D(gaussian->muNum, gaussian->n);
	VectorDouble_destructArray1D(gaussian->muDenom, gaussian->n);
        MatrixDouble_destructArray1D(gaussian->cov, gaussian->n);
	MatrixDouble_destructArray1D(gaussian->covNum, gaussian->n);
	free(gaussian->covDenom);
	free(gaussian->weights);
	free(gaussian->weightNum);
        free(gaussian);
}

/** 
 * Constructs a gaussian object for which the diagonal covariance entries 
 * are initialized to the initial mean values
 *
 * @param mu	An array of double vectors containing the initial mean vectors
 * @param n	The number of mixture components
 * @return gaussian
 */
Gaussian* Gaussian_constructSpecial(VectorDouble** mu, int n){
	// Allocate and initialize a special covariance matrix
	// It would be a diagonal matrix with autovariances same as means
        MatrixDouble** cov = MatrixDouble_constructArray1D(n, mu[0]->dim, mu[0]->dim);
	for (int m=0; m < n; m++){
		for (int j=0; j < mu[m]->dim; j++){
			cov[m]->data[j][j] = mu[m]->data[j];
        	}
	}
	return Gaussian_construct(mu, cov, n);
}


/**
 * Makes a transition matrix with uniform probabilities
 *
 * @param dim	The dimension of the transition matrix supposedly 
 * 		should be equal to the number of states + 1
 * @return trans The transition matrix
 */
MatrixDouble* makeUniformTransition(int dim){
        MatrixDouble* trans = MatrixDouble_construct0(dim, dim);
        for (int i=0; i < dim; i++){
                for (int j=0; j < dim; j++){
                        trans->data[i][j] = 1.0 / dim;
                }
        }
	return trans;
}

/**
 * Makes a transition matrix with probabilites biased toward
 * not changing states
 *
 * @param dim	The dimension of the transition matrix supposedly 
 *              should be equal to the number of states + 1
 * @return trans The transition matrix 
 */
MatrixDouble* makeBiasedTransition(int dim){
        MatrixDouble* trans = MatrixDouble_construct0(dim, dim);
        for (int i=0; i < dim; i++){
                for (int j=0; j < dim; j++){
                        trans->data[i][j] = 1.0 / dim * 0.0001;
                }
		// the diagonal entries have higher probs
		trans->data[i][i] = 1.0 - 0.0001 * (dim - 1.0) / dim;
        }
        return trans;
}


/**
 * Returns the total probabilty of emitting a vector from the given Gaussian
 *
 * @param vec	The emitted vector
 * @param gaussian	The gaussian object
 * @return prob	The emission probability
 */

double getGaussianProb(VectorChar* vec, Gaussian* gaussian){
	double* probs = getGaussianMixtureProbs(vec, gaussian);
	double totProb = 0.0;
	for (int m=0; m < gaussian->n; m++){
		totProb += probs[m];
	}
	free(probs);
	return totProb;
}


/**
 * Returns the probabilties of emitting a vector from the mixture components of the given Gaussian
 *
 * @param vec   The emitted vector
 * @param gaussian      The gaussian object
 * @return prob An array of probabilities for all Gaussian components
 */

// TODO: Support higher dimension if needed. This function only supports dimensions less than or equal to 2.
double* getGaussianMixtureProbs(VectorChar* vec, Gaussian* gaussian){
	uint8_t* x = vec->data;
	int dim = vec->dim;
	double w;
	double* mu;
        MatrixDouble* covCopy;
	double** c;
	double* probs = malloc(gaussian->n * sizeof(double));
	memset(probs, 0, gaussian->n * sizeof(double));
	double det; double a0; double a1;
	// iterate over mixture components
	for (int m=0; m < gaussian->n; m++){
		mu = gaussian->mu[m]->data;
		covCopy = MatrixDouble_copy(gaussian->cov[m]);
		c = covCopy->data;
		w = gaussian->weights[m];
		if (dim == 1){
			probs[m] = w / (c[0][0] * sqrt(2 * PI)) * exp(-0.5 * pow((x[0] - mu[0]) / c[0][0],2));
		}
		else if (dim == 2){
			det = c[1][1] * c[0][0] - c[1][0] * c[0][1];
			// TODO: Find a better solution for zero determinant
                	if(det < 1e-4){
				fprintf(stderr, "The determinant is very low! %.3e\n", det);
				c[1][1] += 0.01;
				c[0][0] += 0.01;
				det = c[1][1] * c[0][0] - c[1][0] * c[0][1];
				fprintf(stderr, "The determinant changed! %.3e\n", det);
			}
			if(det < 0){
				fprintf(stderr, "The determinant has a negative value! %.3e\n", det);
				exit(EXIT_FAILURE);
			}
			a0 = (x[0] - mu[0]) * c[1][1] - (x[1] - mu[1]) * c[1][0];
			a1 = (x[1] - mu[1]) * c[0][0] - (x[0] - mu[0]) * c[0][1];
			probs[m] = w / (sqrt(2.0 * PI) * sqrt(det)) * exp(-0.5 / det * ((x[0] - mu[0]) * a0 + (x[1] - mu[1]) * a1));
		}
		MatrixDouble_destruct(covCopy);
		if(probs[m] != probs[m]){
			double u = -0.5 / det * ((x[0] - mu[0]) * a0 + (x[1] - mu[1]) * a1);
			fprintf(stderr, "prob is NAN exp(%.2e) det=%.2e!\n", u, det);
			fprintf(stderr, "%s\n", MatrixDouble_toString(gaussian->cov[m]));
			fprintf(stderr, "%s\n", VectorDouble_toString(gaussian->mu[m]));
			exit(EXIT_FAILURE);
		}
		if(probs[m] < 1e-40){
			probs[m] = 1e-40;
		}
	}
	/*if(prob < 1e-100){
		double u = -0.5 / det * ((x[0] - mu[0]) * a0 + (x[1] - mu[1]) * a1);
		fprintf(stderr, "#%s\n", VectorDouble_toString(gaussian->mu));
		VectorDouble* muVect = gaussian->mu;
		fprintf(stderr, "x=[%d %d], mu = [%.2f %.2f], det=%.2e, u=%.2e\n", x[0], x[1], muVect->data[0], muVect->data[1], det, u);
	}*/
	return probs;
}

double getExpProb(double x, double lambda){
        double prob = x < 0 ? 0 : lambda * exp(-1 * lambda * x);
        return prob;
}

/**
 * Constructs an HMM object 
 *
 * @param nClasses	The number of region classes like HSats-1,2,and 3
 * @param nComps	The number of components (states) for each region class
 * @param nEmit		The dimension of each emitted vector
 * @param nMixtures	An array that contains the number of mixture components for each Gaussian
 * 			nMixrures[HMM components]
 * @param mu		mu[classes][HMM components][mixture components][emission dim]
 * 			mu[r][c] ---> [mixture components][emission dim]
 * @param hmm		The HMM object
 */

HMM* HMM_construct(int nClasses, int nComps, int nEmit, int* nMixtures, VectorDouble**** mu){
	HMM* model = malloc(sizeof(HMM));
	model->nClasses = nClasses;
	model->nComps = nComps;
	model->nEmit = nEmit;
	model->nMixtures = nMixtures;
	// Constructing the emission Gaussians and set their parameters
	// emit[r][c] is pointing to the Gaussian of the c-th component of the r-th class
	model->emit = (Gaussian***) malloc(nClasses * sizeof(Gaussian**));
	for(int r=0; r < nClasses; r++){
		model->emit[r] = (Gaussian**) malloc(nComps * sizeof(Gaussian*));
		// Initialize the emission parameters of each state with the given vector
		for(int c=0; c < nComps; c++){
			model->emit[r][c] = Gaussian_constructSpecial(mu[r][c], nMixtures[c]);
		}
	}

	// The last row of the transition matrix is pointing to the probabilities of
	// starting the sequence from each possible state
	// The last column of the transition matrix is pointing to the probabilities of
	// ending the sequence with each possible state
	// The dimension of the transition matrix is (#nComps + 1) x (#nComps + 1)
	// Each class of region will have its own transition matrix
	model->trans = (MatrixDouble**) malloc(nClasses * sizeof(MatrixDouble*));
	for(int r=0; r < nClasses; r++){
                model->trans[r] = makeBiasedTransition(nComps + 1);
        }

	// Initialize sufficient stats for estimating transition probs
	model->transCounts = (MatrixDouble**) malloc(nClasses * sizeof(MatrixDouble*));
	for(int r = 0; r < nClasses; r++){
		model->transCounts[r] = MatrixDouble_construct0(nComps + 1, nComps + 1);
	}
	

	return model;
}

/**
 * Destructs an HMM object
 *
 * @param model	an HMM object
 */
void HMM_destruct(HMM* model){

	int nClasses = model->nClasses;
        int nComps = model->nComps;
        int nEmit = model->nEmit;

	//free emit gaussians
        for(int r=0; r < nClasses; r++){
                for(int c=0; c < nComps; c++){
                        Gaussian_destruct(model->emit[r][c]);
                }
		free(model->emit[r]);
        }
	free(model->emit);
	
	//free transition matrices
	for(int r=0; r < nClasses; r++){
                MatrixDouble_destruct(model->trans[r]);
        }
	free(model->trans);

	// free sufficient stats for transition matrix
        for(int r = 0; r < nClasses; r++){
                MatrixDouble_destruct(model->transCounts[r]);
        }
	free(model->transCounts);
	free(model);
}

EM* EM_construct(VectorChar** seqEmit, uint8_t* seqClass, int seqLength, HMM* model){
	assert(seqLength > 0);
	EM* em = malloc(sizeof(EM));
	em->nClasses = model->nClasses;
	em->nComps = model->nComps;
	em->nEmit = model->nEmit;
	// Copy the emission sequence
	em->seqEmit = (VectorChar**) malloc(seqLength * sizeof(VectorChar*));
	for(int i=0; i < seqLength; i++){
		em->seqEmit[i] = VectorChar_copy(seqEmit[i]);
	}
	// Copy the classes sequence
	em->seqClass = (uint8_t*) malloc(seqLength * sizeof(uint8_t));
	memcpy(em->seqClass, seqClass, seqLength);
	// Allocate and initialize forward and backward matrices
	em->seqLength = seqLength;
	em->f = (double**) malloc(seqLength * sizeof(double*));
	em->b = (double**) malloc(seqLength * sizeof(double*));
	for(int i=0; i < seqLength; i++){
                em->f[i] = (double*) malloc(model->nComps * sizeof(double));
		em->b[i] = (double*) malloc(model->nComps * sizeof(double));
		for(int c=0; c < model->nComps; c++){
                        em->f[i][c] = 0.0;
			em->b[i][c] = 0.0;
                }
	}
	// Initialize scale to avoid underflow
	em->scales = (double*) malloc(seqLength * sizeof(double));
	for(int i=0; i < seqLength; i++){
		em->scales[i] = 1.0;
	}
	em->px = -1.0;

	return em;
}

void EM_destruct(EM* em){
	int seqLength = em->seqLength;
	for(int i=0; i < seqLength; i++){
                VectorChar_destruct(em->seqEmit[i]);
        }
	free(em->seqEmit);
	free(em->seqClass);

	for(int i=0; i < seqLength; i++){
		free(em->f[i]);
		free(em->b[i]);
	}
	free(em->f);
	free(em->b);

	free(em->scales);
	free(em);
}

void runForward(HMM* model, EM* em){
	// Get em attributes
	VectorChar** seqEmit = em->seqEmit;
	uint8_t* seqClass = em->seqClass;
	int seqLength = em->seqLength;
	// Get model attributes
	int nComps = model->nComps;
	int nClasses = model->nClasses;
	int nEmit = model->nEmit;
	MatrixDouble** trans = model->trans;
	Gaussian*** emit = model->emit;
	double scale = 0.0;
	// Initialize to zero
	for(int i=0; i < seqLength; i++){
                for(int c=0; c < nComps; c++){
                        em->f[i][c] = 0.0;
                }
        }
	double eProb;
	double tProb;
	// Set the 0-th block of the forward matrix
	scale = 0.0;
	for(int c=0; c < nComps; c++){
		// Emission probability
		eProb = getGaussianProb(seqEmit[0], emit[seqClass[0]][c]);
		// Transition probability
		// trans[]->data[nComps][c] holds the prob of starting with c
		tProb = trans[seqClass[0]]->data[nComps][c];
		// Update forward
		em->f[0][c] = eProb * tProb;
		scale += em->f[0][c];
	}
	em->scales[0] = scale;
	// Scale f
	for(int c=0; c < nComps; c++){
		em->f[0][c] /= scale;
	}
	// Fill the remaining rows of the forward matrix
	for(int i=1; i < seqLength; i++){
		scale = 0.0;
		for(int c2=0; c2 < nComps; c2++){
			for(int c1=0; c1 < nComps; c1++){ // Transition from c1 comp to c2 comp
				if (seqClass[i-1] != seqClass[i]) { // if the class has changed
					// Make the transition prob uniform
					tProb = 1.0 / (nComps + 1); 
				}
				else {
					tProb = trans[seqClass[i]]->data[c1][c2];
				}
                		em->f[i][c2] += (em->f[i-1][c1] * tProb);
			}
			eProb = getGaussianProb(seqEmit[i], emit[seqClass[i]][c2]);
			em->f[i][c2] *= eProb;
			scale += em->f[i][c2];
        	}
		if(scale < 1e-100){
			fprintf(stderr, "scale is very low! %.2e\n", scale);
			for(int c2=0; c2 < nComps; c2++){
				eProb = getGaussianProb(seqEmit[i], emit[seqClass[i]][c2]);
				//fprintf(stderr, "vec=%s,\n mu=%s\nc2=%d, eprob=%.2e\n", VectorChar_toString(seqEmit[i]), VectorDouble_toString(emit[seqClass[i]][c2]), c2, eProb);
			}

		}
		// Scale f
		for(int c=0; c < nComps; c++){
			em->f[i][c] /= scale;
		}
		em->scales[i] = scale;
	}
	// Update P(x)
	em->px = 0;
	for(int c=0; c < nComps; c++){
		tProb = trans[seqClass[seqLength - 1]]->data[c][nComps];
		// data[c][nComps] is the probability of ending with state c
		em->px += em->f[seqLength - 1][c] * tProb;
	}
}


// runForward should be run before runBackward because of scales
void runBackward(HMM* model, EM* em){
	// Get em attributes
        VectorChar** seqEmit = em->seqEmit;
        uint8_t* seqClass = em->seqClass;
        int seqLength = em->seqLength;
        // Get model attributes
        int nComps = model->nComps;
        int nClasses = model->nClasses;
        int nEmit = model->nEmit;
        MatrixDouble** trans = model->trans;
        Gaussian*** emit = model->emit;
        // Initialize to zero
        for(int i=0; i < seqLength; i++){
                for(int c=0; c < nComps; c++){
                        em->b[i][c] = 0.0;
                }
        }
        // Set the last row of the backward matrix
        for(int c=0; c < nComps; c++){
                // trans[]->data[c][nComps] holds the probability of ending with c
                em->b[seqLength - 1][c] = trans[seqClass[seqLength - 1]]->data[c][nComps];
		// Scale b
		em->b[seqLength - 1][c] /= em->scales[seqLength - 1];
		if (em->b[seqLength - 1][c] != em->b[seqLength - 1][c]){
                                fprintf(stderr, "b[seqLength -1] is NAN!\n");
				fprintf(stderr, "scale = %.2e\n", em->scales[seqLength - 1]);
				exit(EXIT_FAILURE);
		}
        }
	double tProb;
	double eProb;
        // Fill the remaining rows of the backward matrix
        for(int i = seqLength - 2; 0 <= i; i--){
                for(int c1=0; c1 < nComps; c1++){
                        for(int c2=0; c2 < nComps; c2++){ // Transition from c1 t c2
				if (seqClass[i] != seqClass[i + 1]){// if the class has changed
					// Make the transition prob uniform
					tProb = 1.0 / (nComps + 1); 
				}
				else{
					tProb = trans[seqClass[i]]->data[c1][c2];
				}
				eProb = getGaussianProb(seqEmit[i+1], emit[seqClass[i+1]][c2]);
                                em->b[i][c1] += (em->b[i+1][c2] * tProb * eProb);
                        }
			// Scale b
                	em->b[i][c1] /= em->scales[i];
			if (em->b[i][c1] != em->b[i][c1]){
				fprintf(stderr, "b is NAN!\n");
				for(int c2=0; c2 < nComps; c2++){
					fprintf(stderr, "c2=%d -> b[i+1][c2] = %.2e, scale = %.2e\n", c2, em->b[i+1][c2], em->scales[i]);
				}
				exit(EXIT_FAILURE);
			}
                }
        }
	// Update P(x)
        em->px = 0;
        for(int c=0; c < nComps; c++){
		// starting with c
		tProb = trans[seqClass[0]]->data[nComps][c];
		eProb = getGaussianProb(seqEmit[0], emit[seqClass[0]][c]);
                em->px += em->b[0][c] * eProb * tProb;
		if (em->px != em->px){
                                fprintf(stderr, "px is NAN!\n");
				exit(EXIT_FAILURE);
		}
        }
}

double* getForward(EM* em, int pos){
        double* f = malloc(em->nComps * sizeof(double));
        for(int c=0; c < em->nComps; c++){
                f[c] = em->f[pos][c];
        }
        return f;
}

double* getBackward(EM* em, int pos){
        double* b = malloc(em->nComps * sizeof(double));
        for(int c=0; c < em->nComps; c++){
                b[c] = em->b[pos][c];
        }
        return b;
}

double* getPosterior(EM* em, int pos){
	double* posterior = malloc(em->nComps * sizeof(double));
	for(int c=0; c < em->nComps; c++){
		posterior[c] = em->f[pos][c] * em->b[pos][c] / em->px * em->scales[pos];
	}
	return posterior;
}


void estimateParameters(HMM* model){
	// Get model attributes
        int nComps = model->nComps;
        int nClasses = model->nClasses;
        int nEmit = model->nEmit;
        MatrixDouble** trans = model->trans;
	// A 1D array of count matrices [nClasses]
        MatrixDouble** transCounts = model->transCounts;

	// Update transition probs
	double terminateProb = 0.00001;
	for(int r=0; r < nClasses; r++){
		fprintf(stderr, "transition counts\nr=%d\n",r);
		fprintf(stderr, "%s\n", MatrixDouble_toString(transCounts[r]));
		for(int c1=0; c1 < nComps; c1++){
			double norm = 0;
			for(int c2=0; c2 < nComps; c2++){
				norm += transCounts[r]->data[c1][c2]; 
			}
			if (norm < 10) continue;
			for(int c2=0; c2 < nComps; c2++){
				trans[r]->data[c1][c2] = transCounts[r]->data[c1][c2] / norm * (1.0 - terminateProb);
			}
			trans[r]->data[c1][nComps] = terminateProb;
		}
		for(int c2=0; c2 < nComps; c2++){
                	trans[r]->data[nComps][c2] = 1.0 / nComps * (1.0 - terminateProb);
                }
		trans[r]->data[nComps][nComps] = terminateProb; 
	}

	double ploidy[5] = {0, 0.5, 1, 1, 2};
	double firstEntry = 0.0; double firstEntryWeight=0.0;
	// Update mu
	for(int r=0; r < nClasses; r++){
		// Sum all numerator and denomerator for estimating the first entry of
		// the mean vector
                for(int c=1; c < nComps; c++){
			Gaussian* gaussian = model->emit[r][c];
			for(int m = 0; m < gaussian->n; m++){
				firstEntry += gaussian->muNum[m]->data[0];
				firstEntryWeight += gaussian->muDenom[m]->data[0];
			}
		}
		// Update the mean vectors
		for(int c=1; c < nComps; c++){ // skipping erroneous component
                        Gaussian* gaussian = model->emit[r][c];
			for(int m = 0; m < gaussian->n; m++){
				if (firstEntryWeight > 50e3){
					gaussian->mu[m]->data[0] = firstEntry / firstEntryWeight * ploidy[c];
				}
				//if (gaussian->muDenom[m]->data[1] > 50e3){
				//	gaussian->mu[m]->data[1] = gaussian->muNum[m]->data[1] / gaussian->muDenom[m]->data[1];
				//}
				if (c == 3 || c== 4) {
					gaussian->mu[m]->data[1] = gaussian->mu[m]->data[0];
				}
			}
		}
	}

	// Update cov
	for(int r=0; r < nClasses; r++){
                for(int c=1; c < nComps; c++){
                        Gaussian* gaussian = model->emit[r][c];
			for(int m = 0; m < gaussian->n; m++){
				if (gaussian->covDenom[m] < 50000) {
					fprintf(stderr, "weight is too small [r=%d, c=%d]: %.1e\n", r, c, gaussian->covDenom[m]);
					continue;
				}
                        	MatrixDouble_copyInPlace(gaussian->cov[m], gaussian->covNum[m]);
                        	MatrixDouble_divideByValue(gaussian->cov[m], gaussian->covDenom[m]);
			}
                }
        }


	// Update mixture weights
	double weightDenom;
	for(int r=0; r < nClasses; r++){
                for(int c=0; c < nComps; c++){
                        Gaussian* gaussian = model->emit[r][c];
			weightDenom = 0.0;
                        for(int m = 0; m < gaussian->n; m++){
                                weightDenom += gaussian->weightNum[m];
                        }
			for(int m = 0; m < gaussian->n; m++){
                                gaussian->weights[m] = gaussian->weightNum[m] / weightDenom;
                        }
                }
        }
}


void resetSufficientStats(HMM* model){
	
	for(int r = 0; r < model->nClasses; r++){
		MatrixDouble_setValue(model->transCounts[r], 0);
	}

	for(int r=0; r < model->nClasses; r++){
		for(int c=0; c < model->nComps; c++){
                        Gaussian* gaussian = model->emit[r][c];
			for(int m = 0; m < gaussian->n; m++){
				VectorDouble_setValue(gaussian->muNum[m], 0.0);
				VectorDouble_setValue(gaussian->muDenom[m], 0.0);
				MatrixDouble_setValue(gaussian->covNum[m], 0.0);
				gaussian->covDenom[m] = 0.0;
				gaussian->weightNum[m] = 0.0;
			}
		}
	}
}

void updateSufficientStats(HMM* model, EM* em){
	// Get em attributes
        VectorChar** seqEmit = em->seqEmit;
        uint8_t* seqClass = em->seqClass;
        int seqLength = em->seqLength;
        // Get model attributes
        int nComps = model->nComps;
        int nClasses = model->nClasses;
        int nEmit = model->nEmit;
        MatrixDouble** trans = model->trans;
	// A 1D array of count matrices [nClasses]
        MatrixDouble** transCounts = model->transCounts;

	double tProb;
	double eProb;
	uint8_t r1;
	uint8_t r2;
	uint8_t r;
	// Update transCounts
	for (int i = 0; i < seqLength - 1; i++){
		r1 = seqClass[i];
		r2 = seqClass[i+1];
		for(int c1 = 0; c1 < nComps; c1++){
			for(int c2 = 0; c2 < nComps; c2++)
				if (r1 == r2){
					tProb = trans[r1]->data[c1][c2];
					eProb = getGaussianProb(seqEmit[i+1], model->emit[r2][c2]);
					double u = em->f[i][c1] * em->b[i+1][c2] / em->px * tProb * eProb;
					transCounts[r1]->data[c1][c2] += u;
					if (u != u){
						fprintf(stderr, "Updating transcounts NAN observed\n");
						fprintf(stderr, "f= %2.e\n", em->f[i][c1]);
						fprintf(stderr, "b= %.2e\n", em->b[i+1][c2]);
						fprintf(stderr, "px= %.2e\n", em->px);
						fprintf(stderr, "tprob=%.2e\n", tProb);
						fprintf(stderr, "eprob=%.2e\n", eProb);
						exit(EXIT_FAILURE);
					}
					if (transCounts[r1]->data[c1][c2] != transCounts[r1]->data[c1][c2]){
						fprintf(stderr, "Updating transcounts NAN observed\n");
						exit(EXIT_FAILURE);
					}

				}
		}
	}

	// TODO: Add weightComps[] to make means dependent on each other
	double ploidy[5] = {0, 0.5, 1, 1, 2};
	double w = 0.0;
	double w1 = 0.0;
	double w2 = 0.0;
	// Update the numerators and denomerators for estimating mu
	// Update the numerator for estimating mixture weights 
	for (int i = 0; i < seqLength - 1; i++){
		r = seqClass[i];
		for(int c = 0; c < nComps; c++){
			Gaussian* gaussian = model->emit[r][c];
			double* mixtureProbs = getGaussianMixtureProbs(seqEmit[i], gaussian);
			double totGaussianProb = getGaussianProb(seqEmit[i], gaussian);
			for(int m = 0; m < gaussian->n; m++){
				w1 = mixtureProbs[m] / totGaussianProb;
				w2 = em->f[i][c] * em->b[i][c] / em->px * em->scales[i];
				w = w1 * w2;
				// Update sufficient stats for estimating mean vectors
				// TODO: it is hard-coded here to only two dimensions!
				if (0 < c){// TODO: starting from c=1 means skipping erroneous component
					gaussian->muNum[m]->data[0] +=  w * seqEmit[i]->data[0] / ploidy[c];
					gaussian->muDenom[m]->data[0] += w;
					//gaussian->muNum[m]->data[1] += w * seqEmit[i]->data[1];
                                	//gaussian->muDenom[m]->data[1] += w;
				}
				// Update sufficient stats for estimating mixture weights
				gaussian->weightNum[m] += w;
				// Update sufficient stats for estimating covariance matrices
				for(int j1 = 0; j1 < nEmit; j1++){
                                	for(int j2 = 0; j2 < nEmit; j2++){
                                        	double z1 = seqEmit[i]->data[j1] - gaussian->mu[m]->data[j1];
                                        	double z2 = seqEmit[i]->data[j2] - gaussian->mu[m]->data[j2];
                                        	gaussian->covNum[m]->data[j1][j2] += w * z1 * z2;
                                	}
				}
				gaussian->covDenom[m] += w;
                        }
			free(mixtureProbs);
		}
	}
		
}


Chunk* Chunk_construct1(int chunkLen){
        Chunk* chunk = malloc(sizeof(Chunk));
        chunk->chunkLen = chunkLen;
        chunk->seqEmit = NULL;
        chunk->seqClass = NULL;
        chunk->seqLen = 0;
        chunk->ctg[0] = '\0';
        chunk->ctgLen = 0;
        chunk->s = -1;
        chunk->e = -1;
        chunk->windowLen = -1;
        chunk->windowItr = -1;
        chunk->windowSumEmit = NULL;
        chunk->windowClass = NULL;
        chunk->fileOffset = 0;
        return chunk;
}

Chunk* Chunk_construct3(int chunkLen, int emissionDim, int windowLen){
	Chunk* chunk = malloc(sizeof(Chunk));
	chunk->chunkLen = chunkLen;
	chunk->seqEmit = VectorChar_constructArray1D(chunkLen * 2, emissionDim);
	chunk->seqClass = malloc(chunkLen * 2 * sizeof(uint8_t));
	chunk->seqLen = 0;
	chunk->ctg[0] = '\0';
	chunk->ctgLen = 0;
	chunk->s = -1;
	chunk->e = -1;
	chunk->windowLen = windowLen;
	chunk->windowItr = -1;
	chunk->windowSumEmit = VectorDouble_construct0(emissionDim);
	chunk->windowClass = malloc(windowLen * sizeof(uint8_t));
	chunk->fileOffset = 0;
	return chunk;
}

void Chunk_destruct(Chunk* chunk){
	if(chunk->seqEmit){
		VectorChar_destructArray1D(chunk->seqEmit, 2 * chunk->chunkLen);
	}
	free(chunk->seqClass);
	if(chunk->windowSumEmit){
		VectorDouble_destruct(chunk->windowSumEmit);
	}
        free(chunk->windowClass);
	free(chunk);
}

Batch* Batch_construct(char* covPath, int chunkLen, int nThreads, int emissionDim, int windowLen){
	Batch* batch = malloc(sizeof(Batch));
	batch->nThreads = nThreads;
        batch->chunkLen = chunkLen;
	batch->nEmit = emissionDim;
	batch->threadChunks = (Chunk**) malloc(nThreads * sizeof(Chunk*));
	for(int i = 0; i < nThreads; i++){
		batch->threadChunks[i] = Chunk_construct3(chunkLen, emissionDim, windowLen);
	}
	char covIndexPath[200];
	sprintf(covIndexPath, "%s.index", covPath);
	batch->templateChunks = parseCovIndex(covIndexPath);
	batch->nThreadChunks = 0;
	batch->templateChunkIdx = 0;
	strcpy(batch->covPath, covPath);
	batch->windowLen = windowLen;
	batch->mutex = malloc(sizeof(pthread_mutex_t));
	pthread_mutex_init(batch->mutex, NULL);
	return batch;
}

void Batch_destruct(Batch* batch){
        for(int i = 0; i < batch->nThreads; i++){
                Chunk_destruct(batch->threadChunks[i]);
        }
	free(batch->threadChunks);
	stList_destruct(batch->templateChunks);
	free(batch->mutex);
        free(batch);
}

uint8_t Chunk_getWindowClass(Chunk* chunk){
	uint8_t nClass = maxCharArray(chunk->windowClass, chunk->windowItr + 1);
	int* counts = malloc((nClass + 1) * sizeof(int));
	memset(counts, 0, (nClass + 1) * sizeof(int));
	for(int i = 0; i < chunk->windowItr + 1; i++){
		counts[chunk->windowClass[i]] += 1;
	}
	uint8_t classMax = 0;
	int maxCount = counts[0];
	for(int i = 1; i < nClass + 1; i++){
                if (maxCount < counts[i]){
			classMax = i;
			maxCount = counts[i];
		}
        }
	free(counts);
	return classMax;
}

int Chunk_addWindow(Chunk* chunk){
	if (chunk->windowItr == -1) return 0;
	VectorDouble* windowSum = chunk->windowSumEmit;
	VectorDouble_divideByValue(windowSum, chunk->windowItr + 1); // take average
	for(int j=0; j < windowSum->dim; j++){
        	chunk->seqEmit[chunk->seqLen]->data[j] = round(windowSum->data[j]); // windowSum is actually average here
	}
	VectorDouble_setValue(chunk->windowSumEmit, 0.0);
	chunk->seqClass[chunk->seqLen] = Chunk_getWindowClass(chunk);
        chunk->seqLen += 1;
        chunk->windowItr = -1;
	return 1;
}


int Chunk_addBlock(Chunk* chunk, Block_t* block){
	// check if the contig name matches
	assert(strcmp(chunk->ctg, block->ctg) == 0);
	int nucLen = chunk->windowLen * chunk->seqLen + chunk->windowItr + 1;
	// check if the overlap starts exactly one base after where the already added sequence ends
	//fprintf(stderr, "%d + %d * %d + %d + 1 == %d\n", chunk->s , chunk->windowLen, chunk->seqLen, chunk->windowItr, max(block->s, chunk->s));
	assert(chunk->s + nucLen  == max(block->s, chunk->s));
	int nucLenToAdd = min(block->e, chunk->e) - max(block->s, chunk->s) + 1;
	if (nucLenToAdd <= 0) return 0;
        VectorDouble* windowSum;
	for (int i = 0; i < nucLenToAdd; i++){
		chunk->windowItr += 1;
                chunk->windowItr %= chunk->windowLen;
		// emitted vector at each location
                windowSum = chunk->windowSumEmit;
                for (int j=0; j < block->attrbsLen - 1; j++){
                    // the attrbs in the block include emission entries and class of region
                    windowSum->data[j] += atoi(block->attrbs[j]);
                }
		// the first attribute right after emission is the class of region
                chunk->windowClass[chunk->windowItr] = atoi(block->attrbs[block->attrbsLen - 1]);
		if (chunk->windowItr == chunk->windowLen - 1){ // the window is fully iterated
			Chunk_addWindow(chunk);
                }
         }
	 return nucLenToAdd;
}


void writeCovIndex(stList* chunks, char* indexPath){
	FILE* filePtr = fopen(indexPath, "w+");
        if (filePtr == NULL){
                printf("Couldn't open %s\n", indexPath);
                exit(EXIT_FAILURE);
        }
	assert( 0 < stList_length(chunks));
	Chunk* chunk = stList_get(chunks, 0);
	fprintf(filePtr, "%d\n", chunk->chunkLen);
	for(int i = 0; i < stList_length(chunks); i++){
		chunk = stList_get(chunks, i);
		fprintf(filePtr, "%s\t%d\t%d\t%d\t%ld\n", chunk->ctg, 
				                          chunk->ctgLen, 
					                  chunk->s,
					                  chunk->e,
					                  chunk->fileOffset);
	}
	fflush(filePtr);
	fclose(filePtr);
	fprintf(stderr, "Index file (%s) is written \n", indexPath);
}

stList* createCovIndex(char* covPath, int chunkLen){
	FILE* filePtr = fopen(covPath, "r+");
        if (filePtr == NULL){
                printf("Couldn't open %s\n", covPath);
                exit(EXIT_FAILURE);
        }
        Block_t* block = Block_construct(COV);
	uint64_t preFileOffset = ftell(filePtr);
	Chunk* chunk = NULL; Chunk* preChunk= NULL;
	stList* chunks = stList_construct3(0, (void(*)(void *))Chunk_destruct);
	while(Block_next(filePtr, block) == 1) {
                block->s -= 1; block->e -= 1; // make them 0-based
                if(chunk == NULL || strcmp(chunk->ctg, block->ctg) != 0){ // contig has changed
                        chunk = Chunk_construct1(chunkLen);
                        chunk->s = 0;
                        chunk->e = block->ctgLen < 2 * chunkLen ? block->ctgLen - 1 : chunkLen - 1;
			strcpy(chunk->ctg, block->ctg);
			chunk->ctgLen = block->ctgLen;
			chunk->fileOffset = preFileOffset;
			stList_append(chunks, chunk);
			//fprintf(stderr, "%d\t%s\t%d\t%d\t%d\t%ld\n", stList_length(chunks), chunk->ctg, chunk->ctgLen, chunk->s, chunk->e, chunk->fileOffset);
                }
		else if (chunk->e <= block->e && chunk->e < block->ctgLen - 1){ // has reached the end of the chunk
			preChunk = chunk;
			chunk = Chunk_construct1(chunkLen);
			chunk->s = preChunk->e + 1;
                        chunk->e = block->ctgLen < preChunk->e + 2 * chunkLen ?  block->ctgLen - 1 : preChunk->e + chunkLen;
			strcpy(chunk->ctg, block->ctg);
			chunk->ctgLen = block->ctgLen;
			if(chunk->s <= block->e)
				chunk->fileOffset = preFileOffset;
			else
				chunk->fileOffset = ftell(filePtr);
			stList_append(chunks, chunk);
			//fprintf(stderr, "%d\t%s\t%d\t%d\t%d\t%ld\n", stList_length(chunks), chunk->ctg, chunk->ctgLen, chunk->s, chunk->e, chunk->fileOffset);
		}
		preFileOffset = ftell(filePtr);
        }
	Block_destruct(block);
        return chunks; // The code should reach here when there is no more block left in the cov file
}

stList* parseCovIndex(char* covIndexPath){
        FILE* filePtr = fopen(covIndexPath, "r+");
        if (filePtr == NULL){
                printf("Couldn't open %s\n", covIndexPath);
                exit(EXIT_FAILURE);
        }
	size_t len = 0;
    	char* line = NULL;
    	char* token;
    	ssize_t read;
	stList* chunks = stList_construct3(0, (void(*)(void *))Chunk_destruct);
	read = getline(&line, &len, filePtr);
	int chunkLen = atoi(line);
	while((read = getline(&line, &len, filePtr)) != -1){
		Chunk* chunk = Chunk_construct1(chunkLen);
		// contig name
		token = strtok(line, "\t");
        	strcpy(chunk->ctg, token);
		// contig length
        	token = strtok(NULL, "\t");
        	chunk->ctgLen = atoi(token);
		// start
		token = strtok(NULL, "\t");
		chunk->s = atoi(token);
		// end
		token = strtok(NULL, "\t");
                chunk->e = atoi(token);
		// file offset
		token = strtok(NULL, "\t");
                chunk->fileOffset = atol(token);
		// add chunk to the list
		stList_append(chunks, chunk);
	}
	return chunks;
}


int Batch_readThreadChunks(Batch* batch){
	if(batch->templateChunkIdx == stList_length(batch->templateChunks)){
		return 0;
	}
	batch->nThreadChunks = 0;
	pthread_t* tids = malloc(batch->nThreads * sizeof(pthread_t));
	int nRunningThreads = 0;
	int nRemainingTemplateChunks = stList_length(batch->templateChunks) - batch->templateChunkIdx;
	while(nRunningThreads < min(batch->nThreads, nRemainingTemplateChunks)){
		pthread_create(&tids[nRunningThreads], NULL, Batch_readNextChunk, (void*) batch);
		//Batch_readNextChunk((void*) batch);
		nRunningThreads += 1;
	}
	for(int t=0; t < nRunningThreads; t++){
		assert(pthread_join(tids[t], NULL) == 0 );
	}
	free(tids);
	return 1;
}

// It may return three numbers:
// 1 means that one new chunk is added successfully
// 0 means that the coverage file has no more chunk to add
int Batch_readNextChunk(void* batch_){
	Batch* batch = (Batch*) batch_;
	// Lock the mutex to stop other threads from reading a new chunk
	pthread_mutex_lock(batch->mutex);
	Chunk* templateChunk = stList_get(batch->templateChunks, batch->templateChunkIdx);
	// Construct a block for iteration
	Block_t* block = Block_construct(COV);
	strcpy(block->ctg, templateChunk->ctg);
	block->ctgLen = templateChunk->ctgLen;
	// Open cov file and jump to the first block of the chunk
	FILE* filePtr = fopen(batch->covPath, "r");
        if (filePtr == NULL){
                printf("Couldn't open %s\n", batch->covPath);
                exit(EXIT_FAILURE);
        }
	fseek(filePtr, templateChunk->fileOffset, SEEK_SET);
	// Get the chunk that has to be filled with the emitted sequence
	// and set start and end coordinates and also the contig name from the given template chunk
	Chunk* chunk = batch->threadChunks[batch->nThreadChunks];
	batch->templateChunkIdx += 1;
        batch->nThreadChunks += 1;
	// Unlock the mutex to let other threads fetch a chunk to read
        pthread_mutex_unlock(batch->mutex);
        chunk->seqLen = 0;
        chunk->windowItr = -1;
        chunk->s = templateChunk->s;
        chunk->e = templateChunk->e;
        strcpy(chunk->ctg, templateChunk->ctg);
        chunk->ctgLen = templateChunk->ctgLen;
	// iterate over the blocks in the cov file
	while(Block_next(filePtr, block) == 1) {
		block->s -= 1; block->e -= 1; // make them 0-based
		// if the block overlaps with the chunk
		if (chunk->s <= block->e && block->s <= chunk->e){
			assert(0 < Chunk_addBlock(chunk, block));
		}
		if (chunk->e <= block->e){ // chunk is read completely
			if (chunk->windowItr != -1){ // handle a partially iterated window
				Chunk_addWindow(chunk);
			}
			Block_destruct(block);
			fclose(filePtr);
			//pthread_mutex_unlock(batch->mutex);
			return 1;
		}
	}
	Block_destruct(block);
	fclose(filePtr);
	//pthread_mutex_unlock(batch->mutex);
	return 0; // The code should reach here when there is no more block left in the cov file
}

