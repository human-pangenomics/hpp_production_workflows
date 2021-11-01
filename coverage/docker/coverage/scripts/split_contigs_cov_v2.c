#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include "sonLib.h"



int min(int a, int b){
	return a < b ? a : b;
}

int max(int a, int b){
        return b < a ? a : b;
}

typedef struct{
	char* contigName;
	int contigLength;
	int start; //inclusive 1-based
	int end; //inclusive 1-based
}Segment;


Segment* Segment_construct(char* contigName, int length, int start, int end){
	Segment* sg = malloc(sizeof(Segment));
	sg->contigName = malloc(50);
	strcpy(sg->contigName, contigName);
	sg->contigLength = length;
	sg->start = start;
	sg->end = end;
}

void Segment_destruct(void* segment){
	Segment* sg = segment;
        free(sg->contigName);
	sg->contigName = NULL;
	free(sg);
}



stList* getSegments(char* faiPath, int segmentSize){
	stList* segments = stList_construct3(0, Segment_destruct);
	char contigName[50]; int s, e, contigSize;
	FILE* faif = fopen(faiPath, "r");
	size_t len = 0;
    	char* line = NULL;
    	char* token;
	int n;
	while(getline(&line, &len, faif) != -1) {
    		token = strtok(line, "\t");
        	strcpy(contigName, token);
        	token = strtok(NULL, "\t");
        	contigSize = atoi(token);
		n = contigSize / segmentSize;
		for(int i=0; i < n; i++){
			s = i * segmentSize + 1;
			e = ( s + 2 * segmentSize - 1 <= contigSize) ? (i + 1) * segmentSize : contigSize;
			stList_append(segments, Segment_construct(contigName, contigSize, s, e)); //inclusive 1-based
		}
		if (n == 0){
			s = 1;
			e = contigSize;
			stList_append(segments, Segment_construct(contigName, contigSize, s, e)); //inclusive 1-based
		}
	}
	return segments;
}

int getBlockTypeIndex(float* probArray){
    float max = -1;
    int index = -1;
    for(int i = 2; i < 6; i++){
        if(max < probArray[i]){
	    index = i;
	    max = probArray[i];
	}
    }
    return index - 2;
}

int readNextBlock(FILE* fileReader, char** contig, int*contigLength, int* blockStart, int* blockEnd, int* coverage){
    size_t len = 0;
    char* line = NULL;
    char* token;
    int start, end;
    ssize_t read = getline(&line, &len, fileReader);
    if (read == -1){
         *contig = NULL;
         *blockStart = -1;
	 *blockEnd = -1;
	 *coverage = -1;
	 return 0;
    }
    if (line[0] == '>') {
	token = strtok(line, " ");
        strcpy(*contig, token + 1); // skip '>' and copy
	token = strtok(NULL, " ");
	*contigLength = atoi(token);
        readNextBlock(fileReader, contig, contigLength, blockStart, blockEnd, coverage);
    }
    else {
	token = strtok(line, "\t");
        *blockStart = atoi(token);
	token = strtok(NULL, "\t");
	*blockEnd = atoi(token);
	token = strtok(NULL, "\t");
        *coverage = atoi(token);
    }
    free(line);
    return 1;
}

void splitCov(char* covPath, char* prefix, stList* segments){
    FILE* fp; FILE* fo = NULL;
    char outputPath[200];
    
    fp = fopen(covPath, "r");
    
    if (fp == NULL)
        exit(EXIT_FAILURE);

    char* contig = malloc(50);
    int contigLength = 0;
    char* preContig = malloc(50);
    int blockStart=0;
    int start=0, end=0, cov=0;
    int sg_idx = 0;
    Segment* sg = stList_get(segments, sg_idx);
    printf("%s\t%d\t%d\n", sg->contigName, sg->start, sg->end);
    sprintf(outputPath, "%s.%s_%d_%d.cov", prefix, sg->contigName, sg->start, sg->end);
    fo = fopen(outputPath, "w");
    while (readNextBlock(fp, &contig, &contigLength, &start, &end, &cov) == 1) {
	if ((strcmp(preContig, contig) != 0) ||  sg->end <= end){
		if ((strcmp(preContig, contig) != 0)) fprintf(fo, ">%s %d\n", contig, contigLength);
		//iterate over segments untill a the end of the currect coverage block is covered
		while (sg && strcmp(sg->contigName, contig) == 0 && sg->end <= end)  {
			fprintf(fo, "%d\t%d\t%d\n", max(start, sg->start), sg->end, cov);
			sg_idx++;
			if(sg_idx < stList_length(segments)) {
			       	sg = stList_get(segments, sg_idx);
			        printf("%s\t%d\t%d\n", sg->contigName, sg->start, sg->end);
				fflush(fo);
                        	fclose(fo);
				sprintf(outputPath, "%s.%s_%d_%d.cov", prefix, sg->contigName, sg->start, sg->end);
				fprintf(fo, ">%s %d\n", sg->contigName, sg->contigLength);
                        	fo = fopen(outputPath, "w");
			}
			else {
			       	sg = NULL;
				break;
			}
		}
		if (sg && strcmp(sg->contigName, contig) == 0 && sg->start < end && end < sg->end) { // write the first part of the next segment
			fprintf(fo, "%d\t%d\t%d\n", sg->start, end, cov);
		}
		strcpy(preContig, contig);
	}
	if (sg && strcmp(sg->contigName, contig) == 0 && sg->start < start && end < sg->end)  {
		fprintf(fo, "%d\t%d\t%d\n", start, end, cov);
	}
    }
    fflush(fo);
    fclose(fo);
    fclose(fp);
}

int main(int argc, char *argv[]) {
   int c;
   int segmentSize=15e6;
   char* faiPath;
   char* covPath;
   char* prefix;
   char *program;
   (program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
   while (~(c=getopt(argc, argv, "c:p:f:s:h"))) {
		switch (c) {
			case 'c':
                                covPath = optarg;
                                break;
			case 'p':
                                prefix = optarg;
                                break;
			case 'f':
                                faiPath = optarg;
                                break;
			case 's':
                                segmentSize = atoi(optarg);
                                break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
help:	
				fprintf(stderr, "\nUsage: %s  -f <FAI_FILE> -c <COVERAGE> -p <PREFIX> -s <SEGMENT_SIZE> \n", program);
				fprintf(stderr, "Options:\n");
				fprintf(stderr, "         -c         coverage file\n");
				fprintf(stderr, "         -f         fai file\n");
				fprintf(stderr, "         -p         prefix for the output cov files\n");
				fprintf(stderr, "         -s         segment size[default : 15Mb]\n");
				return 1;	
		}		
   }
   stList* segments = getSegments(faiPath, segmentSize);
   splitCov(covPath, prefix, segments);
   stList_destruct(segments);
   return 0;
}
