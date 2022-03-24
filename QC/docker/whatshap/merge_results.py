import sys
import collections

# column positions
STATS_CHROM_IDX=1
STATS_FILENAME_IDX=2
PAIRWISE_CHROM_IDX=1
PAIRWISE_FILENAME_IDX=4
PAIRWISE_COVERED_VARIANTS_IDX=7
PAIRWISE_ALL_ASSESSED_PAIRS_IDX=8
PAIRWISE_ALL_SWITCHES_IDX=9
PAIRWISE_ALL_SWITCH_RATE_IDX=10
PAIRWISE_BLOCKWISE_HAMMING_IDX=13
PAIRWISE_BLOCKWISE_HAMMING_RATE_IDX=14

# ensure proper useage
assert(len(sys.argv) == 3)

# parse stats file
stats_header = None
stats = collections.defaultdict(lambda: dict())
with open(sys.argv[1]) as fin:
    for line in fin:
        line = line.strip()
        parts = line.split("\t")
        chrom = parts[STATS_CHROM_IDX]
        filename = parts[STATS_FILENAME_IDX]
        if stats_header is None:
            assert(chrom == "chromosome")
            assert(filename == "file_name")
            stats_header = line
            continue
        assert(chrom not in stats[filename])
        stats[filename][chrom] = line

# parse pairwise file
pairwise_header = None
pairwise = collections.defaultdict(lambda: dict())
with open(sys.argv[2]) as fin:
    for line in fin:
        line = line.strip()
        parts = line.split("\t")
        chrom = parts[PAIRWISE_CHROM_IDX]
        filename = parts[PAIRWISE_FILENAME_IDX]
        if pairwise_header is None:
            assert(chrom == "chromosome")
            assert(filename == "file_name0")
            assert(parts[PAIRWISE_COVERED_VARIANTS_IDX] == "covered_variants")
            assert(parts[PAIRWISE_ALL_ASSESSED_PAIRS_IDX] == "all_assessed_pairs")
            assert(parts[PAIRWISE_ALL_SWITCHES_IDX] == "all_switches")
            assert(parts[PAIRWISE_ALL_SWITCH_RATE_IDX] == "all_switch_rate")
            assert(parts[PAIRWISE_BLOCKWISE_HAMMING_IDX] == "blockwise_hamming")
            assert(parts[PAIRWISE_BLOCKWISE_HAMMING_RATE_IDX] == "blockwise_hamming_rate")
            pairwise_header = line
            continue
        assert(chrom not in pairwise[filename])
        pairwise[filename][chrom] = line

# sanity check
assert (stats_header is not None and pairwise_header is not None)

# calcluate pairwise stats for ALL
for file in pairwise.keys():
    covered_variants = sum(map(lambda x: int(x.split("\t")[PAIRWISE_COVERED_VARIANTS_IDX]), pairwise[file].values()))
    all_assessed_pairs = sum(map(lambda x: int(x.split("\t")[PAIRWISE_ALL_ASSESSED_PAIRS_IDX]), pairwise[file].values()))
    all_switches = sum(map(lambda x: int(x.split("\t")[PAIRWISE_ALL_SWITCHES_IDX]), pairwise[file].values()))
    blockwise_hamming = sum(map(lambda x: int(x.split("\t")[PAIRWISE_BLOCKWISE_HAMMING_IDX]), pairwise[file].values()))
    all_switch_rate = all_switches / all_assessed_pairs
    blockwise_hamming_rate = blockwise_hamming / covered_variants
    pairwise_all_calculated_values = {
        PAIRWISE_CHROM_IDX: "ALL",
        PAIRWISE_FILENAME_IDX: file,
        PAIRWISE_COVERED_VARIANTS_IDX: str(covered_variants),
        PAIRWISE_ALL_ASSESSED_PAIRS_IDX: str(all_assessed_pairs),
        PAIRWISE_ALL_SWITCHES_IDX: str(all_switches),
        PAIRWISE_ALL_SWITCH_RATE_IDX: str(all_switch_rate),
        PAIRWISE_BLOCKWISE_HAMMING_IDX: str(blockwise_hamming),
        PAIRWISE_BLOCKWISE_HAMMING_RATE_IDX: str(blockwise_hamming_rate)
    }
    calculated_all_line = "\t".join(["." if x not in pairwise_all_calculated_values else pairwise_all_calculated_values[x] for x in range(len(pairwise_header.split('\t')))])
    pairwise[file]["ALL"] = calculated_all_line

# print
print(stats_header+"\t"+pairwise_header)
for file in sorted(list(stats.keys() & pairwise.keys())):
    for chrom in sorted(list(stats[file].keys() & pairwise[file].keys())):
        print(stats[file][chrom]+"\t"+pairwise[file][chrom])
