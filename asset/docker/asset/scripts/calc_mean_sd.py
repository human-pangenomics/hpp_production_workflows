import sys
import argparse
import math


def main():
    parser = argparse.ArgumentParser(description='Get the coverage file and calculate mean and sd')
    parser.add_argument('--coverageInput', type=str,
                    help='Input coverage file')
    parser.add_argument('--meanOutput', type=str,
                    help='File path to save coverage mean')
    parser.add_argument('--sdOutput', type=str,
                    help='File path to save coverage standard deviation')
    args = parser.parse_args()
    coveragePath = args.coverageInput
    meanOutputPath = args.meanOutput
    sdOutputPath = args.sdOutput

    # Calculate the mean
    totalSize = 0
    totalCoverage = 0
    with open(coveragePath,'r') as f:
        for line in f:
            if line.startswith(">"):
                continue
            attrbs = line.strip().split()
            start = int(attrbs[0])
            end = int(attrbs[1])
            coverage = int(attrbs[2])
            size = (end - start + 1)
            totalSize += size
            totalCoverage += coverage * size
    coverageMean = totalCoverage / totalSize

    # Calculate the standard deviation
    totalSize = 0
    totalSquaredCoverage = 0
    with open(coveragePath,'r') as f:
        for line in f:
            if line.startswith(">"):
                continue
            attrbs = line.strip().split()
            start = int(attrbs[0])
            end = int(attrbs[1])
            coverage = int(attrbs[2])
            size = (end - start + 1)
            # Exclude very large coverages in calculating standard deviation
            if (coverage < (2.5 * coverageMean)):
                totalSize += size
                totalSquaredCoverage += ((coverage - coverageMean) ** 2) * size
    coverageSD = math.sqrt(totalSquaredCoverage/totalSize)

    with open(meanOutputPath,'w') as f:
        f.write("{:.2f}\n".format(coverageMean))
    with open(sdOutputPath,'w') as f:
        f.write("{:.2f}\n".format(coverageSD))




main()
