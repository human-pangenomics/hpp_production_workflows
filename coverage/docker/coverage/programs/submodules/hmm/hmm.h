#include "block_it.h"
#include "data_types.h"
#include <stdio.h>
#include "sonLib.h"

#ifndef HMM_H
#define HMM_H

typedef struct Gaussian{
        VectorDouble** mu; // Mean vector
        MatrixDouble** cov; // Covariance matrix
	double* weights; // weights of components
	int n; // number of mixture components
	VectorDouble** muNum; // numerator for estimating mean
	VectorDouble** muDenom; // denomerator for estimating mean
	MatrixDouble** covNum; // numerator for estimating cov
	double* covDenom; // denomerator for estimating cov
	double* weightNum; // numerator for estimating mixture weights
} Gaussian;


typedef struct HMM{
        // For each class and each component we have a distinct Gaussian modeling the emission
        Gaussian*** emit;
        // The last row of the transition matrix (trans[nComps][]) holds the starting probs
        // The last column of the transition matrix (trans[][nComps]) holds the terminating probs
        // An array of transition matrices. Each transition matrix has the dimension (#nComps + 1) x (#nComps + 1)
        // And the number of matrices should be equal to nClasses
        MatrixDouble** trans;
        // Saving the counts for updating transitions
        // A 1D array of count matrices [nClasses]
        MatrixDouble** transCounts;
	// Dimension attributes
        int nClasses; // Number of classes like non-HSAT, HSAT1, ...
        int nComps; // Number of components like erroneous, haploid, ...
        int nEmit; // Dimension of each emission
	int* nMixtures; // number of mixture components for each HMM states
} HMM;


typedef struct EM{
        VectorChar** seqEmit; // the sequence of emissions; each emission is a vector of length nEmit
        uint8_t* seqClass; // the sequnce of classes
        int seqLength;
	int nClasses; // Number of classes like non-HSAT, HSAT1, ...
        int nComps; // Number of components like erroneous, haploid, ...
        int nEmit; // Dimension of each emission
        double** f; // Forward matrix: #SEQ_LEN x #nComps
        double** b; // Backward matrix: #SEQ_LEN x #nComps
        double px; // P(x)
        double* scales;
} EM;

typedef struct Chunk{
	// 2 * chunkLen is the maximum size for a chunk
	// the last chunk of a contig is most of the times 
	// longer than chunkLen and shorter than 2 * chunkLen
	VectorChar** seqEmit; // [2 * chunkLen] * [nEmit] the emitted sequence
        uint8_t* seqClass; // [2 * chunkLen] the sequence of class of regions
	char ctg[50]; // the name of the contig this chunk is located 
	int ctgLen; // the length of the contig
	int s; // start location of the chunk on the contig 0-based
	int e; // end location of the chunk on the contig 0-based
	int seqLen; // the length of the sequence that has been added to this chunk so far
	int chunkLen; // size of array is 2 * chunkLen
	int windowLen;
	int windowItr;
	VectorDouble* windowSumEmit;
	uint8_t* windowClass;
	uint64_t fileOffset; // the number of offset bytes to reach the first block of this chunk
} Chunk;


typedef struct Batch{
	// Constant attributes
	int nThreads;
	int chunkLen;
	int nEmit; // The dimension of the emitted vector at each location
	int windowLen;
	// Attributes that may change
	char covPath[200];
	Chunk** threadChunks; // There are at most nThreads chunks in the batch
	stList* templateChunks; // The list of all chunks (no seq) in the cov file
	int templateChunkIdx;
	int nThreadChunks; // The number of chunks added to this batch so far
	pthread_mutex_t* mutex;
} Batch;

Gaussian* Gaussian_construct(VectorDouble** mu, MatrixDouble** cov, int nMixtures);

void Gaussian_destruct(Gaussian* gaussian);

Gaussian* Gaussian_constructSpecial(VectorDouble** mu, int nMixtures);

MatrixDouble* makeUniformTransition(int dim);

double getGaussianProb(VectorChar* vec, Gaussian* gaussian);

double* getGaussianMixtureProbs(VectorChar* vec, Gaussian* gaussian);

HMM* HMM_construct(int nClasses, int nComps, int nEmit, int* nMixtures, VectorDouble**** mu);

void HMM_destruct(HMM* model);

EM* EM_construct(VectorChar** seqEmit, uint8_t* seqClass, int seqLength, HMM* model);

void EM_destruct(EM* em);

void runForward(HMM* model, EM* em);

void runBackward(HMM* model, EM* em);

double* getForward(EM* em, int pos);

double* getBackward(EM* em, int pos);

double* getPosterior(EM* em, int pos);

void estimateParameters(HMM* model);

void resetSufficientStats(HMM* model);

void updateSufficientStats(HMM* model, EM* em);

Chunk* Chunk_construct1(int chunkLen);

Chunk* Chunk_construct3(int chunkLen, int emissionDim, int windowLen);

stList* createCovIndex(char* covPath, int chunkLen);

stList* parseCovIndex(char* covIndexPath);

void writeCovIndex(stList* chunks, char* indexPath);

int Chunk_addBlock(Chunk* chunk, Block_t* block);

Batch* Batch_construct(char* covPath, int chunkLen, int nThreads, int emissionDim, int windowLen);

void Batch_destruct(Batch* batch);

int Batch_readThreadChunks(Batch* batch);

int Batch_readNextChunk(void* batch_);

#endif
