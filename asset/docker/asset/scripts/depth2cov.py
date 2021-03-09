import sys
import argparse


def main():
    parser = argparse.ArgumentParser(description='Compress the output of samtools depth -aa by saving continuous blocks with constant base-level coverage')
    parser.add_argument('--depth', type=str,
                    help='3-column output of samtools depth -aa')
    parser.add_argument('--fai', type=str,
                    help='Genome index in .fai format')
    parser.add_argument('--output', type=str,
                    help='Output file in .cov format')
    args = parser.parse_args()
    depthPath = args.depth
    faiPath = args.fai
    outputPath = args.output

    # Make a dictionary for keeping the size of each contig
    contigSizes = {}
    with open(faiPath,'r') as f:
        for line in f:
            attrbs = line.strip().split()
            contigName = attrbs[0]
            contigSize = int(attrbs[1])
            contigSizes[contigName] = contigSize
    with open(depthPath,'r') as f_in:
        with open(outputPath,'w') as f_out:
            # Read the first position and its coverage
            attrbs = f_in.readline().strip().split()
            contigName = attrbs[0]
            startPos = int(attrbs[1])
            pos = startPos
            coverage = int(attrbs[2])
            f_out.write(">{} {}\n".format(contigName, contigSizes[contigName]))
            for line in f_in:
                attrbs = line.strip().split()
                contigNameNew = attrbs[0]
                posNew = int(attrbs[1])
                coverageNew = int(attrbs[2])
                if contigNameNew != contigName: 
                    # If a new contig is starting here, write the previous coverage
                    f_out.write("{}\t{}\t{}\n".format(startPos, pos, coverage))
                    # Write the name of the new contig and its size
                    f_out.write(">{} {}\n".format(contigNameNew, contigSizes[contigNameNew]))
                    pos = posNew
                    startPos = pos
                    contigName = contigNameNew
                    coverage = coverageNew
                    # If contig and coverage is not changed update the position and continue reading
                elif coverage == coverageNew:
                    pos = posNew
                else:
                    # If contig is the same and coverage changed write the previous coverage and corresponding block
                    f_out.write("{}\t{}\t{}\n".format(startPos, pos, coverage))
                    pos = posNew
                    startPos = pos
                    coverage = coverageNew
            # Write the coverage of the last block
            f_out.write("{}\t{}\t{}\n".format(startPos, pos, coverage))
                    


main()
