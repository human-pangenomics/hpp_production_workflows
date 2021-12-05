## Phasing the long reads aligned to the diploid assemblies

### Motivation
When we align long reads (HiFi or ONT at this time!) to the diploid assembly of the same sample, the reads coming from homozygous regions may
not be aligned to the correct haplotype by the aligner. Other possible locations of such reads are usually reported in the secondary alignments. Finding 
the correct haplotype becomes even harder if our assemblies are erroneous. Breaks and indel errors in the assembly may mislead the aligner. 

In the Figure below you can see an example of two haplotypes that are highly similar and only different in one base. One of these haplotypes is assembled 
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
the other alignments of the same read. A marker that is a mismatch in all of the alignments is removed immediately. Among the remaining markers
we perform two main filtering:

#### 1. 

