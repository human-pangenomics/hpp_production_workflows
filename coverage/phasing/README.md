## Phasing the long reads aligned to the diploid assemblies

### Motivation
When we align long reads (HiFi or ONT at this time!) to the diploid assembly of the same sample, the reads coming from homozygous regions may
not be aligned to the correct haplotype by the aligner. Other possible locations of such reads are usually reported in the secondary alignments. Finding 
the correct haplotype becomes even harder if our assemblies are erroneous. Breaks and indel errors in the assembly may mislead the aligner. 

In the Figure below you can see an example of two haplotypes that are highly similar and different only in one base. One of these haplotypes is assembled 
correctly and the other one has a long indel or a break point (could be split into two separate contigs). If we align a read from the first haplotype to
the diploid assembly, the aligner may align the read to the both haplotypes and report the one to the false haplotype as the primary alignment. This happens
because of the error in the assembly. Even if the assembly is perfect it may also happen for highly similar haplotypes since 
the aligner chooses one haplotype randomly.
<img src="https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/phasing/images/phase_reads_1.png" width="700" height="400">

Given a set of secondary alignments (most of the time only one) and one primary alignment we want to make sure if 
the primary one is to the correct haplotype and if it is not find the correct one among the secondary alignments.

### Approach

One way to achieve this aim is to find the single-base markers that can navigate us toward the correct haplotype. An example of such marker is shown
in the top figure by green and red bars. All mismatched bases in the alignments of a single read form the initial set of candidate markers. After 
projecting the intial markers into the coordinate of the read (instead of reference) it is possible to determine the match/mismatch status of a marker in 
the other alignments of the same read. A marker that is a mismatch in all of the alignments is removed immediately. 

Among the remaining markers we perform two main filtering:

#### 1. Filter Markers within insertions

If a marker, which is a mismatch in at least one alignment, appears within an insertion in another alignment we remove it in this step. We remove it 
because it can be either a misassembly on the haplotype that contains the insertion or an error on the read. So that marker can be misleading on either cases. In the figure below you can see an example of such a marker (shown with a blue bar).

<img src="https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/phasing/images/phase_reads_2.png" width="700" height="400">

#### 2. Filter markers with low BAQ

After filtering the markers within insertions the remaining markers are all either match or mismatch in any of the alignments. 
The other issue is that sometimes the alignment around a marker is not reliable due to the errors in the reads or the assembly especially homopolymer run errors. To measure the reliability of the alignment around each marker we calculate 
Base Alignment Quality ( BAQ ). This is an adjustment to the raw base quality of the marker and works by realiging the reads to where they were 
already aligned to. This realignment is performed through a banded HMM which parameters have to be tuned before hand. 
There are three parameters that have to be tuned for each sequencing platform; gap opening probability, gap extenstion probability and the bandwidth of the HMM.

After tuning the parameters based on the platform (which is either HiFi or ONT here) we filtered the markers with BAQs lower than a specific threshold (20 for HiFi and 10 for ONT). 

There are two points worth noting:
- BAQ cannot be larger than the raw base quality reported by the sequencer (or base caller)
- For better performance the bases far from the markers (farther than 500) are not included in the BAQ calculation.

More information about BAQ can be found in this paper. BAQ calculation is already implemented in htslib and that implmentation has been imported to this pipeline.

<img src="https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/phasing/images/phase_reads_3.png" width="700" height="400">

#### Select the best alignment using Marker Consistency Score

After finding the confident markers we calculate marker consistency score for each alignment. For calcultaing this score for each alignment we take the markers that appeared as mismatches on that alignment and take a summation of their BAQ values with a negative sign. After calcualating the marker consistency score for all primary and secondary alignments of the same read we select the one with the largest score (usually zero is the largest) and report it as the alignment to the correct haplotype. If the selected alignment is primary we do nothing. You can see an example of calculating marker consistency score in the figure above. 

Two heuristics are applied for increasing specificity:
- If the selected alignment is having a score lower than a specific threshold (in this pipeline `-50`) we do not report it as the best one.
- If the selected alignment is secondary it should be at least `20` (`10` for ONT) units larger than the score of the primary alignment otherwise it is not reported as the correct alignment.
