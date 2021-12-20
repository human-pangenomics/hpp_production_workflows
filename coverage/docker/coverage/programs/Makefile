CC=gcc
INC=-I /home/apps/sonLib/C/inc/ -I /usr/local/include/htslib -I submodules/cigar_it -I submodules/common
STATIC_LIBS=/home/apps/sonLib/lib/sonLib.a
CCFLAGS=$(STATIC_LIBS) $(INC) -lm -lhts
BIN_DIR=bin
SRC_DIR=src

CIGAR_IT_OBJ=submodules/cigar_it/cigar_it.o
COMMON_OBJ=submodules/common/common.o
SRCS=$(wildcard ${SRC_DIR}/*.c)
BINS=$(patsubst ${SRC_DIR}/%.c,${BIN_DIR}/%,$(SRCS))

.PHONY: all

all: ${COMMON_OBJ} ${CIGAR_IT_OBJ} ${BIN_DIR} $(BINS) ${BIN_DIR}/phase_reads_debug

bin/%: src/%.c
	$(CC) -o $@ $< $(CIGAR_IT_OBJ) $(COMMON_OBJ) $(CCFLAGS)

bin/phase_reads_debug: src/phase_reads.c
	$(CC) -o$@ $< $(CIGAR_IT_OBJ) $(COMMON_OBJ) -DDEBUG $(CCFLAGS)

${BIN_DIR}:
	@mkdir -p bin

$(CIGAR_IT_OBJ): submodules/cigar_it/cigar_it.c submodules/cigar_it/cigar_it.h $(COMMON_OBJ)
	$(CC) -c $< -o $@ $(COMMON_OBJ) $(CCFLAGS)

$(COMMON_OBJ): submodules/common/common.c submodules/common/common.h
	$(CC) -c $< -o $@ $(CCFLAGS)

clean:
	/bin/rm -rf $(BINS)
	/bin/rm -rf $(CIGAR_IT_OBJ)
	/bin/rm -rf $(COMMON_OBJ)
