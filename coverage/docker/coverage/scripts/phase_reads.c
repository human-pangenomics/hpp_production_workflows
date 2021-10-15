#include <getopt.h>
#include "sam.h"
#include "faidx.h"
#include <time.h>
#include "bgzf.h"
#include <regex.h>
#include "sonLib.h"
#include <assert.h>
#include <math.h>
#include <float.h>

#define ARRAY_SIZE(arr) (sizeof((arr)) / sizeof((arr)[0]))

#define DEBUG

#ifdef DEBUG
//void DEBUG_PRINT(const char *, ...);
void DEBUG_PRINT(const char *fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
}

#else
static inline void DEBUG_PRINT(const char *fmt, ...) {};
#endif


/*! @typedef
 @abstract Structure for a marker information
 @field alignment_idx	The index of the alignment (could be secondary or primary) from where this marker is taken
 @field read_pos_f	The position of the marker in the read (All positions should be computed w.r.t to
			the forward strand for the inconsistency between negative and positive alignments)
 @field base_q	The base quality of the marker (could be the raw quality or any processed form e.g BAQ)
 @field is_match		True if the marker base is equal to the reference and False otherwise
 */
typedef struct {
	int32_t alignment_idx;
	int32_t read_pos_f;
	int32_t base_idx;
	int32_t base_q;
	bool is_match;
}ptMarker;


/// Create a ptMarker structure
/**
   @param alignment_idx   The index of the alignment (could be secondary or primary) from where this marker is taken
   @param read_pos_f        The position of the marker in the read (All positions should be computed w.r.t to
   			  the forward strand for the inconsistency between negative and positive alignments)
   @param base_q    The base quality of the marker (could be the raw quality or any processed form e.g BAQ)
   @param is_match         True if the marker base is equal to the reference and False otherwise
   @return Pointer to a new marker struct on success; NULL on failure
   The ptMarker struct returned by a successful call should be freed
   via free() when it is no longer needed.
 */
ptMarker* ptMarker_construct(int32_t alignment_idx, int32_t base_idx, int32_t read_pos_f, uint8_t base_q, bool is_match){
	ptMarker* marker = (ptMarker*) malloc(sizeof(ptMarker));
	marker->alignment_idx = alignment_idx;
        marker->is_match = is_match;
        marker->read_pos_f = read_pos_f;
	marker->base_idx = base_idx;
        marker->base_q = base_q;
	return marker;
}

ptMarker* ptMarker_copy(ptMarker* src){
        ptMarker* dest = (ptMarker*) malloc(sizeof(ptMarker));
        dest->alignment_idx = src->alignment_idx;
        dest->is_match = src->is_match;
        dest->read_pos_f = src->read_pos_f;
        dest->base_idx = src->base_idx;
        dest->base_q = src->base_q;
        return dest;
}

/*! @typedef
 @abstract Structure for an alignment record and its likelihood
 @field record	Pointer to a bam1_t struct that contains the alignment record
 @field p	The likelihood of being the correct alignment
 */
typedef struct {
        bam1_t* record; 
	double qv; // The probability of being the correct alignment in log scale (-10*log10(P))
	stList* conf_blocks; // confident blocks
}ptAlignment;


typedef struct {
        int rfs; // ref start
	int rfe; // ref end
	int sqs; // seq start
	int sqe; // seq end
}ptBlock;

ptBlock* ptBlock_construct(int rfs, int rfe, int sqs, int sqe){
	ptBlock* block = malloc(sizeof(ptBlock));
	block->rfs = rfs;
	block->rfe = rfe;
	block->sqs = sqs;
	block->sqe = sqe;
	return block;
}

typedef struct {
	// the constant attributes
        uint32_t* cigar;
	int n; // number of operations
	bool is_rev;
	// the attributes below may change in each iteration
	uint8_t op; // current cigar operation
	int len; // length of the current operation
        int idx;
        int rds_f; // read start forward
	int rde_f; // read end forward
	int sqs; // seq start
	int sqe; // seq end
	int rfs; // ref start
	int rfe; // ref end
}ptCigarIt;

ptCigarIt* ptCigarIt_construct(bam1_t* b){
	ptCigarIt* cigar_it = (ptCigarIt*) malloc(sizeof(ptCigarIt));
	cigar_it->cigar = bam_get_cigar(b);
	cigar_it->len = 0;
	cigar_it->op = -1;
	cigar_it->is_rev = bam_is_rev(b);
	cigar_it->n = b->core.n_cigar;
	cigar_it->idx = -1;
	cigar_it->sqs = 0;
	cigar_it->sqe = -1;
	cigar_it->rfs = b->core.pos;
	cigar_it->rfe = b->core.pos - 1;
	int lclip = 0;
	int rclip = 0;
	if (bam_cigar_op(cigar_it->cigar[0]) == BAM_CHARD_CLIP) {
		lclip = bam_cigar_oplen(cigar_it->cigar[0]);
        }
        if (bam_cigar_op(cigar_it->cigar[cigar_it->n - 1]) == BAM_CHARD_CLIP) {
        	rclip = bam_cigar_oplen(cigar_it->cigar[cigar_it->n - 1]);
        }
	cigar_it->rds_f = cigar_it->is_rev ? rclip + lclip + b->core.l_qseq : 0;
	cigar_it->rde_f = cigar_it->is_rev ? rclip + lclip + b->core.l_qseq  - 1 : -1;
	return cigar_it;
}

int ptCigarIt_next(ptCigarIt* cigar_it){
	if (cigar_it->idx == cigar_it->n - 1) return 0;
	uint32_t* cigar = cigar_it->cigar;
	cigar_it->idx += 1;
	int rd_step;
	int sq_step;
	int rf_step;
	cigar_it->op = bam_cigar_op(cigar[cigar_it->idx]);
	cigar_it->len = bam_cigar_oplen(cigar[cigar_it->idx]);
	switch(cigar_it->op) {
		case BAM_CMATCH:
		case BAM_CEQUAL:
		case BAM_CDIFF:
			rd_step = cigar_it->len;
			sq_step = cigar_it->len;
			rf_step = cigar_it->len;	
      			break;
		case BAM_CINS:
			rd_step = cigar_it->len;
                        sq_step = cigar_it->len;
                        rf_step = 0;
                        break;
   		case BAM_CDEL:
			rd_step = 0;
			sq_step = 0;
			rf_step = cigar_it->len;
			break;
		case BAM_CSOFT_CLIP:
			rd_step = cigar_it->len;
			sq_step = cigar_it->len;
			rf_step = 0;
			break;
		case BAM_CHARD_CLIP:
			rd_step = cigar_it->len;
                        sq_step = 0;
                        rf_step = 0;
      			break;
	}
	//shift read interval on forward strand
	if (cigar_it->is_rev) {
		cigar_it->rde_f = cigar_it->rds_f - 1;
		cigar_it->rds_f -= rd_step;
	}
	else {
		cigar_it->rds_f = cigar_it->rde_f + 1;
		cigar_it->rde_f += rd_step;
	}
	//shift seq interval
	cigar_it->sqs = cigar_it->sqe + 1;
	cigar_it->sqe += sq_step;
	//shift ref interval
	cigar_it->rfs = cigar_it->rfe + 1;
	cigar_it->rfe += rf_step;
	return cigar_it->len;
}

ptAlignment* ptAlignment_construct(bam1_t* record, double qv){
	ptAlignment* alignment = (ptAlignment*) malloc(sizeof(ptAlignment));
	bam1_t* record_cpy = bam_init1();
	assert(bam_copy1(record_cpy, record) != NULL);
	alignment->record = record_cpy;
	alignment->qv = qv;
	alignment->conf_blocks = NULL;
	return alignment;
}

void ptAlignment_destruct(ptAlignment* alignment){
	if (alignment->record){
		bam_destroy1(alignment->record);
        	alignment->record = NULL;
	}
	if (alignment->conf_blocks){
		stList_destruct(alignment->conf_blocks);
		alignment->conf_blocks = NULL;
	}
	free(alignment);
}


/// Compare two ptMarker structures
/**
   @param a	Pointer to the first ptMarker structure
   @param b	Pointer to the second ptMarker structure
   @return 	Nagative when a should be before b and positive otherwise
   @note	The comparison is first based on their positions in the read and if both given markers
		are located at the same position then we look at the alignment_idx 
 */
int ptMarker_cmp(const void *a, const void *b){
        ptMarker* m1 = (ptMarker*) a;
        ptMarker* m2 = (ptMarker*) b;
	if (m1->read_pos_f == m2->read_pos_f) {
		return m1->alignment_idx - m2->alignment_idx;
	}
	else {
		return m1->read_pos_f - m2->read_pos_f;
	}
}

ptMarker* ptMarker_construct_match(ptAlignment** alignments, int32_t alignment_idx, int32_t read_pos_f){
	bam1_t* record = alignments[alignment_idx]->record;
        uint8_t* q = bam_get_qual(record);
	ptMarker* match;
	uint32_t* cigar = bam_get_cigar(record);
	int lclip=0;
	int rclip=0;
	if (bam_cigar_op(cigar[0]) == BAM_CHARD_CLIP) {
                lclip = bam_cigar_oplen(cigar[0]);
        }
        if (bam_cigar_op(cigar[record->core.n_cigar - 1]) == BAM_CHARD_CLIP) {
                rclip = bam_cigar_oplen(cigar[record->core.n_cigar - 1]);
        }
	// construct the match marker
	if (bam_is_rev(record)){
		match = ptMarker_construct(alignment_idx,
					   record->core.l_qseq + rclip - read_pos_f - 1,
                                           read_pos_f,
                                           q[record->core.l_qseq + rclip - read_pos_f - 1],
                                           true);
        }
        else {
		match = ptMarker_construct(alignment_idx,
				           read_pos_f - lclip,
                                           read_pos_f,
                                           q[read_pos_f - lclip],
                                           true);
        }
	return match;
}


stList* filter_lowq_markers(stList* markers, int threshold){
	int markers_len = stList_length(markers);
        stList* keep_markers = stList_construct3(0, free);
	ptMarker* pre_marker =NULL;
	ptMarker* marker;
	bool keep_flag = true;
	int idx_s = 0;
	for(int j=0; j < markers_len; j++){
		marker = stList_get(markers, j);
		if (pre_marker && marker->read_pos_f != pre_marker->read_pos_f){
			if (keep_flag){
				for(int k=idx_s; k<j; k++){
					// copy all markers with the same read_pos_f and add them
					stList_append(keep_markers, ptMarker_copy(stList_get(markers, k)));
				}
			}
			idx_s = j;
			keep_flag = true;
		}
		if(marker->base_q < threshold) keep_flag=false;
		pre_marker = marker;
	}
	if (keep_flag){
		for(int k=idx_s; k<markers_len; k++){
			stList_append(keep_markers, ptMarker_copy(stList_get(markers, k)));
                }
        }
	return keep_markers;
}

stList* filter_ins_markers(stList* markers, ptAlignment** alignments, int alignments_len){
	int markers_len = stList_length(markers);
	bool* markers_keep_flags = (bool*) malloc(markers_len * sizeof(bool));
	for(int i=0; i< markers_len; i++) markers_keep_flags[i] =true;
	bam1_t* b;
	ptMarker* marker;
	ptCigarIt* cigar_it;
	int j; // marker index
        for(int i=0; i < alignments_len; i++){
		b = alignments[i]->record;
		j = bam_is_rev(b) ? stList_length(markers) - 1 : 0;
		marker = stList_get(markers, j);
		cigar_it = ptCigarIt_construct(b);
		//iterate over cigars
		while(ptCigarIt_next(cigar_it)){
			//check if marker is located within the cigar operation
			while(marker &&
			      (cigar_it->rds_f <= marker->read_pos_f) && 
			      (cigar_it->rde_f >= marker->read_pos_f)){
				//check if the marker is located within an insertion
				if(cigar_it->op == BAM_CINS ||
				   cigar_it->op == BAM_CSOFT_CLIP ||
				   cigar_it->op == BAM_CHARD_CLIP){
					markers_keep_flags[j] = false;
				}
				j += bam_is_rev(b) ? -1 : 1;
				marker = ((j < markers_len) && (j >= 0)) ? stList_get(markers, j) : NULL;
			}	
		}//end of iteration over cigar ops
	}// end of iteration over alignments
	// make a list for saving remaining markers
	stList* keep_markers = stList_construct3(0, free);
	for(j=0; j < markers_len; j++){
		if(markers_keep_flags[j]){
			marker = stList_get(markers, j);
			stList_append(keep_markers, ptMarker_copy(marker));
		}
	}
	return keep_markers;
}

void sort_and_fill_markers(stList* markers, ptAlignment** alignments, int alignments_len){
	stList_sort(markers, ptMarker_cmp);
	int idx=0;
	ptMarker* match;
	ptMarker* marker;
	ptMarker* pre_marker = NULL;
	int markers_len = stList_length(markers);
	for(int64_t i=0; i < markers_len; i++){
		marker = stList_get(markers, i);
		assert(!marker->is_match);
		if (pre_marker  &&
		    pre_marker->read_pos_f < marker->read_pos_f){
			for (int64_t j=idx; j < alignments_len; j++){
                        	match = ptMarker_construct_match(alignments, j, pre_marker->read_pos_f);
                        	stList_append(markers, match);
                	}
			idx = 0;
		}
		pre_marker = marker;
		for (int64_t j=idx; j < marker->alignment_idx; j++){
			match = ptMarker_construct_match(alignments, j, marker->read_pos_f);
                	stList_append(markers, match);
        	}
		idx = marker->alignment_idx + 1;
	}
	// add the last matches if there still exist to be added
	for (int64_t j=idx; j < alignments_len; j++){
		match = ptMarker_construct_match(alignments, j, marker->read_pos_f);
        	stList_append(markers, match);
        }
	stList_sort(markers, ptMarker_cmp); // sort again to locate the matches in their correct positions
}

double revert_quality(uint8_t q){
	if (q >= 93) return 0;
	if (q == 0) return 93;
	double p = 1 - pow(10, (double)q / -10);
	double rev_q = -10 * log(p);
	return rev_q;
}

void calc_likelihood(stList* markers, ptAlignment** alignments){
        ptMarker* marker;
        for(int64_t j=0; j < stList_length(markers); j++){
                marker =  stList_get(markers, j);
		if (marker->is_match){
			alignments[marker->alignment_idx]->qv -= revert_quality(marker->base_q);
		}
		else{
			alignments[marker->alignment_idx]->qv -= marker->base_q ;
		}
	}
}

int get_best_record(ptAlignment** alignments, int alignments_len, double min_qv, double prim_margin){
	assert(alignments_len > 0);
	if (alignments_len == 1) return 0;
	double max_qv = -DBL_MAX;
	int max_idx = -1;
	double prim_qv;
	int prim_idx;
	for(int i=0; i < alignments_len; i++){
		if ((alignments[i]->record->core.flag & BAM_FSECONDARY) == 0){
			prim_idx = i;
			prim_qv = alignments[i]->qv;
		}
		else if (max_qv < alignments[i]->qv){
			max_idx = i;
			max_qv = alignments[i]->qv;
		}
	}
	if (max_qv < prim_qv + prim_margin || max_qv < min_qv) return prim_idx;
	else return max_idx;
}

stList* find_confident_blocks(bam1_t* b, int threshold){
	stList* conf_blocks = stList_construct3(0, free);
	ptCigarIt* cigar_it = ptCigarIt_construct(b);
	int conf_sqs = 0; // the seq start of the confident block
	int conf_rfs = b->core.pos; // the ref start of the confident block
	ptBlock* block = NULL;
	while(ptCigarIt_next(cigar_it)){
		switch(cigar_it->op) {
                	case BAM_CINS:
			case BAM_CDEL:
                        	if(cigar_it->len > threshold && 
				   conf_sqs < cigar_it->sqs &&
				   conf_rfs < cigar_it->rfs){
					block = ptBlock_construct(conf_rfs, cigar_it->rfs - 1,
								  conf_sqs, cigar_it->sqs - 1);
					stList_append(conf_blocks, block);
				}
				if(cigar_it->len > threshold){
					conf_sqs = cigar_it->sqe + 1;
                                	conf_rfs = cigar_it->rfe + 1;
				}
                        	break;
                	case BAM_CSOFT_CLIP:
			case BAM_CHARD_CLIP:
                        	if(conf_sqs < cigar_it->sqs &&
				   conf_rfs < cigar_it->rfs){
                                        block = ptBlock_construct(conf_rfs, cigar_it->rfs - 1,
                                                                  conf_sqs, cigar_it->sqs - 1);
					stList_append(conf_blocks, block);
				}
                                conf_sqs = cigar_it->sqe + 1;
                                conf_rfs = cigar_it->rfe + 1;
                        	break;
        	}

	}
	// when the last block is not terminated with SOFT or HARD clipping
	if(conf_sqs <= cigar_it->sqe){
		//printf("LAST (%d-th) BLOCK: %d\t%d\n", stList_length(confident_blocks) - 1, confident_seq_start, cigar_it->seq_end);
		block = ptBlock_construct(conf_rfs, cigar_it->rfe,
                                          conf_sqs, cigar_it->sqe);
		stList_append(conf_blocks, block);
	}
	return conf_blocks;
}


void set_confident_blocks(ptAlignment** alignments, int alignments_len, int threshold){
	for(int i=0; i < alignments_len; i++){
		alignments[i]->conf_blocks = find_confident_blocks(alignments[i]->record, threshold);
	}
}


void calc_local_baq(const faidx_t* fai, const char* contig_name, ptAlignment* alignment, int alignment_idx, stList* markers, double conf_d, double conf_e, double conf_bw){

	uint8_t* tseq; // translated seq A=>0,C=>1,G=>2,T=>3,other=>4
	uint8_t* tref; // translated ref
	uint8_t* bq;
	uint8_t* block_qual;
	int* state;
	uint8_t* q; // Probability of incorrect alignment from probaln_glocal()
	probaln_par_t conf = { conf_d, conf_e, conf_bw };
	char* ref;
	char* reg;
	bam1_t* b = alignment->record;
	uint8_t* seq = bam_get_seq(b);
	uint8_t* qual = bam_get_qual(b);
	stList* blocks = alignment->conf_blocks;
	int markers_len = stList_length(markers);
	int j = bam_is_rev(b) ? stList_length(markers) - 1 : 0;
	ptMarker* marker = stList_get(markers, j);
	ptCigarIt* cigar_it = ptCigarIt_construct(b);
	ptBlock* block;
	int ref_len; int seq_len;
	DEBUG_PRINT("\t\t\t### Number of Blocks: %ld\n", stList_length(blocks));
	for(int i=0; i < stList_length(blocks); i++){
		block = stList_get(blocks, i);
		DEBUG_PRINT("\t\t\t### block#%d: seq[%d:%d] ref[%s:%d-%d]\n", i, block->sqs, block->sqe, contig_name, block->rfs, block->rfe);
		while((cigar_it->sqs < block->sqs) || (cigar_it->rfs < block->rfs)){
			if(ptCigarIt_next(cigar_it) == 0) break;
		}// end of while cigar_it->seq_start should be now equal to block->seq_start
		assert(cigar_it->sqs == block->sqs);
		assert(cigar_it->rfs == block->rfs);
		while(marker &&
                      ((block->sqs > marker->base_idx) ||
		        (marker->alignment_idx != alignment_idx))){
			j += bam_is_rev(alignment->record) ? -1 : 1;
                        marker = ((j < markers_len) && (j >= 0)) ? stList_get(markers, j) : NULL;
		}// end of while the marker is now the first one located within the current block (or after)
		if(marker &&
                   (block->sqe >= marker->base_idx)){
			seq_len = block->sqe - block->sqs + 1;
			ref_len = block->rfe - block->rfs + 1;
			DEBUG_PRINT("\t\t\t### seq_len:%d\tref_len:%d\n", seq_len, ref_len);
			if (seq_len < 5 || ref_len < 5) continue; // skip small blocks
			//allocate read sequence
			tseq = (uint8_t*) malloc(seq_len);
			for(int k=0; k < seq_len; k++)
				tseq[k] = seq_nt16_int[bam_seqi(seq, block->sqs + k)];
			//allocate reference sequenc
			tref = (uint8_t*) malloc(ref_len);
			reg = malloc(200);
			memset(reg,'\0',200);
			sprintf(reg, "{%s}:%d-%d", contig_name, block->rfs + 1 , block->rfe + 1);
			int len;
			ref = fai_fetch(fai, reg, &len);
			assert(len == (block->rfe - block->rfs + 1));
			for(int k=0; k < ref_len; k++)
                                tref[k] = seq_nt16_int[seq_nt16_table[(unsigned char)ref[k]]];
			//allocate the quality of this confident block
			block_qual = (uint8_t*) malloc(seq_len);
                        for(int k=0; k < seq_len; k++)
				block_qual[k] = qual[block->sqs + k];
			//allocate neccessary arrays for HMM BAQ
			state = (int*) malloc((block->sqe - block->sqs + 1) * sizeof(int));
			q = (uint8_t*) malloc(block->sqe - block->sqs + 1);
			//DEBUG_PRINT("Starting local BAQ :))\n");	
			if (probaln_glocal(tref, block->rfe - block->rfs + 1, 
					   tseq, block->sqe - block->sqs + 1,
					   block_qual, &conf, state, q) == INT_MIN) {
            			fprintf(stderr, "probaln_glocal ERROR\n");
				fprintf(stderr, "%s:%d-%d", contig_name, block->rfs + 1 , block->rfe + 1);

			}
			//DEBUG_PRINT("local BAQ Finished:))\n");
			// the state and q are now updated if there is any marker located within the block
			//apply BAQ to the quality array
			bq = (uint8_t*) malloc(seq_len);
			memcpy(bq, block_qual, seq_len);
			int x;
			int y;
			while(cigar_it->sqe <= block->sqe && cigar_it->rfe <= block->rfe){
				x = cigar_it->rfs - block->rfs;
				y = cigar_it->sqs - block->sqs;
				//DEBUG_PRINT("rf:%d\t%d\n", cigar_it->rfe, block->rfe);
				//DEBUG_PRINT("sq:%d\t%d\n", cigar_it->sqe, block->sqe);
				if (cigar_it->op == BAM_CMATCH || 
			    	    cigar_it->op == BAM_CEQUAL || 
			            cigar_it->op == BAM_CDIFF) {
					for (int t = y; t < (y + cigar_it->len); t++) {
						assert(t < seq_len);
                        			if (((state[t]&3) != 0) || (state[t]>>2 != x + (t - y))) bq[t] = 0;
						else bq[t] = bq[t] < q[t] ? bq[t] : q[t];
                    			}
				}
				if(ptCigarIt_next(cigar_it) == 0) break;
			}
		
			for(int k=0; k < seq_len; k++) qual[block->sqs + k] = bq[k] < 94 ? bq[k] : 93;
			
			free(tseq);
			free(tref);
			free(state);
			free(q);
			free(bq);
			free(block_qual);
			free(ref);
			free(reg);
		}
	}
}

void calc_update_baq_all(const faidx_t* fai, 
		        ptAlignment** alignments, int alignments_len, 
			stList* markers, 
			const sam_hdr_t *h,
			double conf_d, double conf_e, double conf_b){
	const char* contig_name;
	int tid;
	for(int i=0; i < alignments_len; i++){
		tid = alignments[i]->record->core.tid;
		contig_name = sam_hdr_tid2name(h, tid);
		calc_local_baq(fai, contig_name, alignments[i], i, markers, conf_d, conf_e, conf_b);
	}
	ptMarker* marker;
	uint8_t* q;
	//DEBUG_PRINT("\t\t## start updating marker base qualities\n");
	for(int i=0; i < stList_length(markers); i++){
		marker = stList_get(markers, i);
		q = bam_get_qual(alignments[marker->alignment_idx]->record);
		marker->base_q = q[marker->base_idx];
	}
}

bool contain_supp(ptAlignment** alignments, int alignments_len){
	for(int i=0;i < alignments_len; i++){
		if (alignments[i]->record->core.flag & BAM_FSUPPLEMENTARY) return 1;
	}
	return 0;
}
int main(int argc, char *argv[]){
	int c;
	bool baq_flag=false;
	double conf_d;
	double conf_e;
	double conf_b;
	char* inputPath;
	char* outputPath;
	char* fastaPath;
   	char *program;
   	(program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
   	while (~(c=getopt(argc, argv, "i:f:o:q:d:e:b:h"))) {
		switch (c) {
                        case 'i':
                                inputPath = optarg;
                                break;
                        case 'o':
                                outputPath = optarg;
                                break;
			case 'f':
				fastaPath = optarg;
				break;
			case 'q':
				baq_flag = true;
				break;
			case 'd':
                                conf_d = atof(optarg);
                                break;
			case 'e':
                                conf_e = atof(optarg);
                                break;
			case 'b':
                                conf_b = atof(optarg);
                                break;
                        default:
                                if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
			help:
                                fprintf(stderr, "\nUsage: %s  -i <INPUT_BAM> -f <FASTA> -o <output_BAM> \n", program);
                                fprintf(stderr, "Options:\n");
                                fprintf(stderr, "         -i         input bam file\n");
				fprintf(stderr, "         -f         input fasta file\n");
                                fprintf(stderr, "         -o         output bam file\n");
				fprintf(stderr, "         -q         optional (use BAQ as base quality)\n");
				fprintf(stderr, "         -d         gap prob\n");
				fprintf(stderr, "         -e         gap extension\n");
				fprintf(stderr, "         -b         dp bandwidth\n");
                                return 1;
		}
	}

	faidx_t* fai = fai_load(fastaPath);
	samFile *fo = sam_open(outputPath, "w");
	samFile* fp = sam_open(inputPath, "r");
	sam_hdr_t* sam_hdr = sam_hdr_read(fp);
	sam_hdr_write(fo, sam_hdr);
	bam1_t* b = bam_init1();
	char read_name[100];
	char read_name_new[100];
	memset(read_name, '\0', 100);
	memset(read_name_new, '\0', 100);
        uint8_t* quality;
	ptMarker* marker;
	stList* markers = stList_construct3(0, free);
	stList* no_ins_markers = NULL;
	stList* highq_markers = NULL;
	int alignments_len=0;
	int32_t read_pos;
	ptAlignment** alignments = (ptAlignment**) malloc(100 * sizeof(ptAlignment*)); //Assuming that no more than 100 alignments we have per read
	int bytes_read;
	const char* contig_name;
	bool supp_flag = false;
	while(true) {
		bytes_read = sam_read1(fp, sam_hdr, b);
		if (bytes_read > - 1){
			strcpy(read_name_new, bam_get_qname(b));
			if (read_name[0] == '\0') {
				strcpy(read_name, read_name_new);
			}
		}
		// If read name has changed or file is finished
		if ((strcmp(read_name_new, read_name) != 0) || (bytes_read <= -1)){
			// If we have at least one marker and at least two alignments
			// then we can decide which one is the best alignment
			if ((stList_length(markers) > 0) && (alignments_len > 1)){
				DEBUG_PRINT("@@ READ NAME: %s\n\t$ Number of alignments: %d\t Read l_qseq: %d\n", read_name, alignments_len, alignments[0]->record->core.l_qseq);
				DEBUG_PRINT("\t# Set confident blocks");
				set_confident_blocks(alignments, alignments_len, 2);
				DEBUG_PRINT("\t$ length of markers: %ld\n\t# Start filtering markers\n", stList_length(markers));
				no_ins_markers = filter_ins_markers(markers, alignments, alignments_len);
				DEBUG_PRINT("\t$ length of confident markers: %ld\n", stList_length(no_ins_markers));
				if(stList_length(no_ins_markers) > 0){
					DEBUG_PRINT("\t# There are more than 0 markers So calc baq and update quality\n");
					sort_and_fill_markers(no_ins_markers, alignments, alignments_len);
					/*for(int t =0 ; t<stList_length(confident_markers); t++) {
                                                ptMarker* m = stList_get(confident_markers, t);
                                                DEBUG_PRINT("MARKER@%d\t%d\t%d\t%d\t%d\t%d\n", t, m->alignment_idx, m->base_idx, m->read_pos_f, m->base_q, m->is_match);
                                        }*/

					if(baq_flag){
						calc_update_baq_all(fai, 
						                    alignments, alignments_len,
						                    no_ins_markers, sam_hdr,
								    conf_d, conf_e, conf_b);
					}
					highq_markers = filter_lowq_markers(no_ins_markers, 20);
					DEBUG_PRINT("@MARKERS:\n");
					for(int t =0 ; t<stList_length(highq_markers); t++) {
                                                ptMarker* m = stList_get(highq_markers, t);
                                                DEBUG_PRINT("MARKER@%d\t%d\t%d\t%d\t%d\t%d\n", t, m->alignment_idx, m->base_idx, m->read_pos_f, m->base_q, m->is_match);
                                        }
					//DEBUG_PRINT("\t# calc likelihood\n");
					calc_likelihood(highq_markers, alignments);
					stList_destruct(highq_markers);
                                	highq_markers = NULL;
				}
				stList_destruct(no_ins_markers);
                                no_ins_markers = NULL;
			}
			if (alignments_len > 0){ // maybe the previous alignment was unmapped
				// get the best alignment and write
				int best_idx = get_best_record(alignments, alignments_len, -100.0, 20.0);
				bam1_t* best = alignments[best_idx]->record;
				// write all chimeric alignments without any change
				if(contain_supp(alignments, alignments_len)){
					for(int i=0; i<alignments_len; i++){
						if (sam_write1(fo, sam_hdr, alignments[i]->record) == -1) {
                                                	fprintf(stderr, "Couldn't write %s\n", read_name);
                                             	}
					}
				}
				else if ((best->core.flag & BAM_FSECONDARY)){
					//make it primary
					printf("$ RECORDS AND QVs for %s:\n", read_name);
					for(int i=0; i<alignments_len; i++){
						if ((alignments[i]->record->core.flag & BAM_FSECONDARY) == 0){
							printf("*\t");
							alignments[i]->record->core.flag |= BAM_FSECONDARY;
							if (sam_write1(fo, sam_hdr, alignments[i]->record) == -1) {
                                        			fprintf(stderr, "Couldn't write %s\n", read_name);
                  					}
						}
						else if (i == best_idx) printf("@\t");
						else printf("!\t");
						contig_name = sam_hdr_tid2name(sam_hdr, alignments[i]->record->core.tid);
						printf("%.2f\t%s\t%ld\n", alignments[i]->qv, contig_name, alignments[i]->record->core.pos);

					}
					printf("\n");
					best->core.flag &= ~BAM_FSECONDARY;
                                        if (sam_write1(fo, sam_hdr, best) == -1) {
                                                fprintf(stderr, "Couldn't write %s\n", read_name);
                                        }
				}
			}
			stList_destruct(markers);
			markers = stList_construct3(0, free);
			for(int i = 0; i < alignments_len; i++){
				ptAlignment_destruct(alignments[i]);
                		alignments[i] = NULL;
        		}
			// initialize for new alignments
			alignments_len = 0;
			strcpy(read_name, read_name_new);
		}
		if (bytes_read <= -1) break; // file is finished so break
		if(b->core.flag & BAM_FUNMAP) continue; // unmapped
		alignments[alignments_len] = ptAlignment_construct(b, 0.0);
		ptCigarIt* cigar_it = ptCigarIt_construct(b);
		quality = bam_get_qual(b);
		while(ptCigarIt_next(cigar_it)){
			if (cigar_it->op == BAM_CDIFF) {
				for(int j=0; j < cigar_it->len; j++){
                                        if (bam_is_rev(b)) {
						marker = ptMarker_construct(alignments_len,
                                                                            cigar_it->sqs + j,
                                                                            cigar_it->rde_f - j,
                                                                            quality[cigar_it->sqs + j],
                                                                            false);
					}
					else{
						marker = ptMarker_construct(alignments_len,
                                                                            cigar_it->sqs + j,
                                                                            cigar_it->rds_f + j,
                                                                            quality[cigar_it->sqs + j],
                                                                            false);
					}
					stList_append(markers, marker);
				}
			}
		}
		alignments_len += 1;
	}
	printf("Done!\n");
	// free memory
	fai_destroy(fai);
	sam_close(fp);
	sam_close(fo);
	bam_destroy1(b);
        sam_hdr_destroy(sam_hdr);
}

//main();
