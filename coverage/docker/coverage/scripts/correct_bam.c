#include <getopt.h>
#include "sam.h"
#include <assert.h>
#include "faidx.h"
#include <time.h>
#include "bgzf.h"
#include <regex.h>
#include "sonLib.h"

typedef struct {
	char*	contig_name;
        int32_t start;
	uint8_t mapq;
}Location;

Location* location_construct(char* contig_name, int start, uint8_t mapq){
	Location* loc = malloc(sizeof(Location));
	loc->contig_name = contig_name;
	loc->start = start;
	loc->mapq = mapq;
	return loc;
}

void location_destruct(void* location){
	Location* loc = location;
	free(loc->contig_name);
	free(loc);
}

stHash* get_phased_read_table(char* phased_reads_path){
	stHash* phased_read_table = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL, location_destruct);
	if(phased_read_table == NULL) return phased_read_table;
	FILE* fp = fopen(phased_reads_path, "r");
	size_t read;
	size_t len;
	char* line = NULL;
	char* read_name;
	char* contig_name;
	int start;
	char* token;
	Location* loc;
	while((read = getline(&line, &len, fp)) != -1){
		line[strlen(line)-1] = '\0';
		if(line[0] == '$'){
			token = strtok(line, "\t");
			token = strtok(NULL, "\t");
			read_name = malloc(strlen(token) + 1);
			strcpy(read_name, token);
		}
		else if (line[0] == '@'){
			token = strtok(line, "\t");
                        token = strtok(NULL, "\t");
			token = strtok(NULL, "\t");
			contig_name = malloc(strlen(token) + 1);
			strcpy(contig_name, token);
			token = strtok(NULL, "\t");
			start = atoi(token); // 0-based
			loc = location_construct(contig_name, start, -1);
			stHash_insert(phased_read_table, read_name, loc);
		}
	}
	return phased_read_table;
}

bool is_prim(stHash* phased_read_table, bam1_t* b, sam_hdr_t* sam_hdr){
	int start = b->core.pos; // 0-based
	const char* contig_name = sam_hdr_tid2name(sam_hdr, b->core.tid);
	char* read_name = bam_get_qname(b);
	Location* loc = stHash_search(phased_read_table, read_name);
	if(loc != NULL){
		if(strcmp(loc->contig_name, contig_name) == 0 && 
		   loc->start == start){
			return true;
			printf("%s\t%s\t%d\n", read_name, loc->contig_name, loc->start);
		}
		else{
			return false;
		}
	}
	return (b->core.flag & BAM_FSECONDARY) == 0; 
}

stHash* get_mapq_table(char* mapq_table_path){
	stHash* mapq_table = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL, (void (*)(void *)) stList_destruct);
	if (mapq_table_path == NULL) return mapq_table;
        FILE* fp = fopen(mapq_table_path, "r");
        size_t read;
        size_t len;
        char* line = NULL;
        char* read_name;
        char* contig_name;
        int start;
	int mapq;
        char* token;
        Location* loc;
	stList* locs;
        while((read = getline(&line, &len, fp)) != -1){
                line[strlen(line)-1] = '\0';
                token = strtok(line, "\t");
                read_name = malloc(strlen(token) + 1);
               	strcpy(read_name, token);
                token = strtok(NULL, "\t");
                contig_name = malloc(strlen(token) + 1);
                strcpy(contig_name, token);
                token = strtok(NULL, "\t");
                start = atoi(token); // 0-based
		token = strtok(NULL, "\t");
                mapq = atoi(token);
                loc = location_construct(contig_name, start, mapq);
		locs = stHash_search(mapq_table, read_name);
		if (locs == NULL){
			locs = stList_construct3(0, location_destruct);
			stHash_insert(mapq_table, read_name, locs);
		}
		stList_append(locs, loc);
        }
        return mapq_table;
}

uint8_t get_mapq(stHash* mapq_table, bam1_t* b, sam_hdr_t* sam_hdr){
        int start = b->core.pos; // 0-based
        const char* contig_name = sam_hdr_tid2name(sam_hdr, b->core.tid);
        char* read_name = bam_get_qname(b);
       	stList* locs = stHash_search(mapq_table, read_name);
	Location* loc;
        if(locs != NULL){
		for(int i = 0; i < stList_length(locs); i++){
			loc = stList_get(locs, i);
			if(strcmp(loc->contig_name, contig_name) == 0 &&
			   loc->start == start){
				return loc->mapq;
			}
		}
        }
        return b->core.qual;
}



int main(int argc, char *argv[]){
        int c;
        char* input_path;
	char* output_path;
	char* phasing_log_path;
	char* mapq_table_path;
	bool primary_only = false;
	bool no_tag = false;
        char *program;
        (program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
        while (~(c=getopt(argc, argv, "i:o:p:m:t:r:h"))) {
                switch (c) {
                        case 'i':
                                input_path = optarg;
                                break;
                        case 'o':
                                output_path = optarg;
                                break;
			case 'p':
                                phasing_log_path = optarg;
                                break;
			case 'm':
				mapq_table_path = optarg;
				break;
			case 't':
				no_tag = true;
				break;
			case 'r':
				primary_only = true;
				break;
                        default:
                                if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
                        help:
                                fprintf(stderr, "\nUsage: %s  -i <INPUT_BAM> -o <OUTPUT_BAM> -p <PHASING_LOG> -m <MAPQ_TABLE>\n", program);
                                fprintf(stderr, "Options:\n");
                                fprintf(stderr, "         -i         input bam file\n");
                                fprintf(stderr, "         -o         output bam file\n");
				fprintf(stderr, "         -p         the phasing log path\n");
				fprintf(stderr, "         -m         the adjusted mapq table path\n");
                                fprintf(stderr, "         -t         output no optional fields\n");
                                fprintf(stderr, "         -m         output only primary alignments\n");
				return 1;
                }
        }
	//open input bam file and read the header
        samFile* fp = sam_open(input_path, "r");
        sam_hdr_t* sam_hdr = sam_hdr_read(fp);
	//open output bam file and write the header
	samFile* fo = sam_open(output_path, "wb");
	sam_hdr_write(fo, sam_hdr);
	//get the table of phased reads (the locations of their primary alignments are stored as values)
	stHash* phased_read_table = get_phased_read_table(phasing_log_path);
	//get the table of reads with adjusted mapq (the locations (+ corresponding mapqs) are stored as values)
	stHashIterator* it = stHash_getIterator(phased_read_table);
	char* read_name;
	Location* loc;
	stHash* mapq_table = get_mapq_table(mapq_table_path);
	//modify mapq and sec/pri tags based on given logs and tables
	bam1_t* b = bam_init1();
	while(sam_read1(fp, sam_hdr, b) > -1) {
		if (is_prim(phased_read_table, b, sam_hdr)){
			b->core.flag &= ~BAM_FSECONDARY; // make it primary
		}
		else{
			if(primary_only) continue;
			b->core.flag |= BAM_FSECONDARY; // make it secondary
		}
		b->core.qual = get_mapq(mapq_table, b, sam_hdr);
		if(no_tag) b->l_data -= bam_get_l_aux(b);
		if (sam_write1(fo, sam_hdr, b) == -1) {
			fprintf(stderr, "Couldn't write %s\n", bam_get_qname(b));
                }
	}
	sam_close(fo);
	sam_close(fp);
        bam_destroy1(b);
        sam_hdr_destroy(sam_hdr);
	stHash_destruct(mapq_table);
	stHash_destruct(phased_read_table);
}
