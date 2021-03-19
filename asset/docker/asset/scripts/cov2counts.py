#from coverage_frequency_distribution import CoverageDistribution
from collections import defaultdict
import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description='Counts the coverage frequencies')
    parser.add_argument('--cov', type=str, help='Coverage file')
    parser.add_argument('--output', type=str, help='Output file in .counts format')
    args = parser.parse_args()
    covPath = args.cov
    outputPath = args.output

    counts = defaultdict(int)
    with open(covPath,"r") as f:
        for line in f:
            if line.startswith(">"):
                print(line.strip())
                continue
            attrbs = line.strip().split()
            start = int(attrbs[0])
            end = int(attrbs[1])
            cov = int(attrbs[2])
            if cov in counts:
                counts[cov] += end - start + 1
            else:
                counts[cov] = end - start + 1

    with open(outputPath,"w") as f:
        for cov in range(len(counts)):
            f.write("{}\t{}\n".format(cov, counts[cov]))


    #cov_model = CoverageDistribution.fit(counts, tol = 1e-4, max_iters=1000)
    #for cov in (10, 20, 30, 40, 50, 60, 70, 80, 90, 100):
    #    prob_err = cov_model.probability_erroneous(cov)
    #   prob_correct = cov_model.probability_haploid(cov)
    #    prob_collapsed = cov_model.probability_collapsed(cov)

    #    print("k-mer frequency {}: prob error {:.4f}, hap {:.4f}, coll {:.4f},".format(cov, prob_err, prob_correct, prob_collapsed))

if __name__ == "__main__":
    main()

