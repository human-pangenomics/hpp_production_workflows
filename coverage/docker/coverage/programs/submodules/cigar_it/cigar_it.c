#include "sam.h"
#include "cigar_it.h"
#include <regex.h>
#include <getopt.h>


void ptCigarIt_destruct(ptCigarIt* cigar_it){
	regfree(&cigar_it->preg_cs);
	free(cigar_it);
}

ptCigarIt* ptCigarIt_construct(bam1_t* b, bool use_cs){
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
	cigar_it->use_cs = use_cs;
	int error;
	if(use_cs){
		cigar_it->offset_cs = 0;
		cigar_it->cs = bam_aux_get(b, "cs") + 1; // skip "Z" at the beginning
		if ((error = regcomp(&cigar_it->preg_cs, CS_PATTERN, REG_EXTENDED)) != 0) {
                	printf("regcomp() failed, returning nonzero (%d)\n", error);
                	exit(1);
        	}
		//printf("%s\n",cigar_it->cs);
	}
	return cigar_it;
}


int ptCigarIt_next_cs(ptCigarIt* cigar_it){
        regmatch_t pm;
	char* cs_shifted = cigar_it->cs + cigar_it->offset_cs;
	int error = regexec(&cigar_it->preg_cs, cs_shifted, 1, &pm, REG_NOTBOL);
	if (error != 0) return 0;
	int rd_step;
        int sq_step;
        int rf_step;
	char match_len[20];
	/*for(int i=pm.rm_so; i<pm.rm_eo; i++){
		printf("%c", cs_shifted[i]);
	}
	printf("\n");*/
        switch(cs_shifted[pm.rm_so]) {
		case ':':
			cigar_it->op = BAM_CEQUAL;
			memcpy(match_len, cs_shifted + 1, pm.rm_eo - pm.rm_so - 1);
			match_len[pm.rm_eo - pm.rm_so - 1] = '\0';
			cigar_it->len = atoi(match_len);
			rd_step = cigar_it->len;
                        sq_step = cigar_it->len;
                        rf_step = cigar_it->len;
                        break;
		case '*':
			cigar_it->op = BAM_CDIFF;
			cigar_it->len = 1;
                        rd_step = cigar_it->len;
                        sq_step = cigar_it->len;
                        rf_step = cigar_it->len;
                        break;
                case '+':
			cigar_it->op = BAM_CINS;
			cigar_it->len = pm.rm_eo - pm.rm_so - 1;
                        rd_step = cigar_it->len;
                        sq_step = cigar_it->len;
                        rf_step = 0;
                        break;
                case '-':
			cigar_it->op = BAM_CDEL;
			cigar_it->len = pm.rm_eo - pm.rm_so - 1;
                        rd_step = 0;
                        sq_step = 0;
                        rf_step = cigar_it->len;
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
	//update offset for cs tag
	cigar_it->offset_cs += pm.rm_eo;
	return cigar_it->len;
}

int ptCigarIt_next(ptCigarIt* cigar_it){
	if (cigar_it->idx == cigar_it->n - 1) return 0;
	uint32_t* cigar = cigar_it->cigar;
	if(cigar_it->use_cs){ // if it is set to use CS tag instead of cigar itself
		bool first_clip = cigar_it->idx == -1 && 
				  (bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP ||
                           	   bam_cigar_op(cigar[0]) == BAM_CHARD_CLIP);
		bool second_clip = cigar_it->idx == 0 &&
                                  (bam_cigar_op(cigar[1]) == BAM_CSOFT_CLIP ||
                                   bam_cigar_op(cigar[1]) == BAM_CHARD_CLIP);
		if(!(first_clip || second_clip)){//if no soft- or hard-clip at the beginning
			int cs_return = ptCigarIt_next_cs(cigar_it);
			// if CS is finished and there exist soft- or hard- clipping 
			// in one operation before the last one
			if(cs_return == 0 &&
			   (bam_cigar_op(cigar[cigar_it->n-2]) == BAM_CSOFT_CLIP ||
                           bam_cigar_op(cigar[cigar_it->n-2]) == BAM_CHARD_CLIP)){
				cigar_it->idx = cigar_it->idx = cigar_it->n-3;
				cigar_it->use_cs = false;
				return ptCigarIt_next(cigar_it);
                        }
			// if CS is finished and there exist soft- or hard- clipping
                        // in the last operation
			else if(cs_return == 0 &&
				(bam_cigar_op(cigar[cigar_it->n-1]) == BAM_CSOFT_CLIP ||
                           	bam_cigar_op(cigar[cigar_it->n-1]) == BAM_CHARD_CLIP )){
                                cigar_it->idx = cigar_it->idx = cigar_it->n-2;
                                cigar_it->use_cs = false;
                                return ptCigarIt_next(cigar_it);
                        }
			else return cs_return;
		}
	}
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
