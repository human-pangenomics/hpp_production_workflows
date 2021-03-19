import sys
import argparse
from collections import defaultdict
import re

def handleCigar(cigarString):
    cigarEvents = re.compile("[0-9]*").split(cigarString)[1:]
    cigarSizes = [int(size) for size in re.compile("H|S|M|I|D").split(cigarString)[:-1]]
    cigarList = [ (event,size) for event, size in zip(cigarEvents,cigarSizes)]
    return cigarList

def findMappedBlocks(cigarList, blocks, chromStart):
    refCurrentPos = chromStart
    cigarStartIdx = 0
    blockIdx = 0
    mappedBlocks = []
    qBlocks = []
    mappedStartPos = None
    contigCurrentPos = 0
    if cigarList[0][0] == 'H' or cigarList[0][0] == 'S':
        #print("#",cigarList[0][0],cigarList[0][1])
        contigCurrentPos = cigarList[0][1]
        cigarStartIdx = 1
        #print("#",blockIdx, blocks)
        while (blockIdx < len(blocks)) and (contigCurrentPos > blocks[blockIdx][1]):
            blockIdx += 1
        if blockIdx >= len(blocks):
            return mappedBlocks, qBlocks
        if contigCurrentPos >= blocks[blockIdx][0]:
            mappedStartPos = refCurrentPos
            qStartPos = contigCurrentPos
    #print(contigCurrentPos,refCurrentPos)
    # iterate over cigar elements and find the blocks
    for cigarEvent, cigarSize in cigarList[cigarStartIdx:]:
        #print("#",cigarEvent,cigarSize)
        # If mismatch or match
        if cigarEvent == 'M':
            contigNextPos = contigCurrentPos + cigarSize
            while (blockIdx < len(blocks)) and (contigNextPos >= blocks[blockIdx][1]):
                if contigCurrentPos < blocks[blockIdx][0]:
                    mappedStartPos = refCurrentPos + blocks[blockIdx][0] - contigCurrentPos
                    qStartPos = blocks[blockIdx][0]
                mappedEndPos = refCurrentPos + blocks[blockIdx][1] - contigCurrentPos
                qEndPos = blocks[blockIdx][1]
                mappedBlocks.append((mappedStartPos, mappedEndPos))
                qBlocks.append((qStartPos, qEndPos))
                blockIdx += 1
            if blockIdx >= len(blocks):
                break
            if (contigCurrentPos < blocks[blockIdx][0]) and (contigNextPos >= blocks[blockIdx][0]):
                mappedStartPos = refCurrentPos + blocks[blockIdx][0] - contigCurrentPos
                qStartPos = blocks[blockIdx][0]
            refCurrentPos += cigarSize
            contigCurrentPos += cigarSize 
        # If insertion
        elif cigarEvent == 'I':
            contigNextPos = contigCurrentPos + cigarSize
            while (blockIdx < len(blocks)) and (contigNextPos >= blocks[blockIdx][1]):
                # If a block is completely within an insertion from start to end ignore it
                if contigCurrentPos < blocks[blockIdx][0]:
                    blockIdx += 1
                    continue
                mappedEndPos = refCurrentPos
                qEndPos = blocks[blockIdx][1]
                mappedBlocks.append((mappedStartPos, refCurrentPos))
                qBlocks.append((qStartPos, qEndPos))
                blockIdx += 1
            if blockIdx >= len(blocks):
                break
            if (contigCurrentPos < blocks[blockIdx][0]) and (contigNextPos >= blocks[blockIdx][0]):
                mappedStartPos = refCurrentPos
                qStartPos = blocks[blockIdx][0]
            contigCurrentPos += cigarSize
        # If deletion
        elif cigarEvent == 'D':
            refCurrentPos += cigarSize
        elif (cigarEvent == 'H') or (cigarEvent == 'S'):
            if (contigCurrentPos >= blocks[blockIdx][0]) and (mappedStartPos < refCurrentPos):
                mappedEndPos = refCurrentPos
                mappedBlocks.append((mappedStartPos, refCurrentPos))
        #print(contigCurrentPos,refCurrentPos)
    return mappedBlocks, qBlocks
                


def main():
    parser = argparse.ArgumentParser(description='Extract blocks')
    parser.add_argument('--sam', type=str,
                    help='sam file')
    parser.add_argument('--bed', type=str,
                    help='bed file')
    parser.add_argument('--output', type=str,
                    help='Output file in bed format')
    args = parser.parse_args()
    samPath = args.sam
    bedPath = args.bed
    outputPath = args.output

    blocks = defaultdict(list)
    with open(bedPath,"r") as f:
        for line in f:
            attrbs = line.strip().split()
            contigName = attrbs[0]
            # start is 0-based in bed format
            start = int(attrbs[1]) + 1
            end = int(attrbs[2])
            blocks[contigName].append((start, end))

    refBlocksAll = []
    qBlocksAll = []
    with open(samPath,"r") as f:
        for line in f:
            attrbs = line.strip().split()
            contigName = attrbs[0]
            chrom = attrbs[2]
            chromStart = int(attrbs[3])
            cigarList = handleCigar(attrbs[5])
            if len(blocks[contigName]) == 0:
                continue
            refBlocks, qBlocks = findMappedBlocks(cigarList, blocks[contigName], chromStart - 1)
            refBlocksAll.extend([(chrom, refBlock[0], refBlock[1]) for refBlock in refBlocks])
            qBlocksAll.extend([(contigName, qBlock[0], qBlock[1]) for qBlock in qBlocks])
    
    with open(outputPath, "w") as f:
        for refBlock in refBlocksAll:
            f.write("{}\t{}\t{}\n".format(refBlock[0], refBlock[1]-1, refBlock[2]))
    #with open("test.bed", "w") as f:
    #    for qBlock in qBlocksAll:
    #        f.write("{}\t{}\t{}\n".format(qBlock[0], qBlock[1], qBlock[2]))




main()
