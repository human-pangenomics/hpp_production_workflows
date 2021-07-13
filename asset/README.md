## Coverage Analysis With Asset For HPRC_PLUS Y1 Assemblies
** This page is written specifically for the asset results provided for the T2T cenSat paper. 


### Overview
This page explains the steps toward flagging the unreliable blocks in HPRC_PLUS Y1 assemblies. Here is a summary of this analysis 
- Align the HiFi reads to each assembly
- Calculate the read coverage of each base in the assembly
- Set the lower- and higher-bound coverage thresholds 
- Run asset to flag the unreliable blocks

### 1. Read Alignment
The HiFi reads are aligned to each assembly with winnowmap v2.0.
Here are the main commands used for producing the alignments (adopted from the [winnowmap docs](https://github.com/marbl/Winnowmap)):
```` 
  # making the k-mer table with meryl
  meryl count k=15 output merylDB asm.fa
  meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt
   # alignment with winnowmap
  winnowmap -W repetitive_k15.txt -ax map-pb asm.fa reads.fq.gz | samtools view -hb > read_alignment.bam
````

### 2. Calculating Depth of Coverage
With `samtools depth` we could calculate the depth of coverage for each base of an assembly. The output of `samtools depth -aa` is like below. (`-aa` 
option allows reporting the bases with no coverage)
````
contig_1  1 0
contig_1  2 1
contig_1  3 1
contig_1  4 1
contig_1  5 2
contig_1  6 2
````
In the order they appear, the columns are showing the contig name, the base coordinate and the coverage. 
Each base has a separate line even if consecutive bases are having the same coverage. 
In order to make this output more efficient we coverted it to the format below
````
>contig_1 6
1 1 0
2 4 1
5 6 2
````
Where the name of each contig appears only once before the first block of that contig and the consecutive bases with the same coverage take only one line.
The number that comes after the name of the contig is the contig size. We added the suffix `.cov` to the files with this format.
Here are the main commands for producing the coverage files:
````
samtools depth -aa -Q 20 alignment.bam > alignment.depth
./depth2cov -d alignment.depth -f asm.fa.fai -o read_alignment.cov
````
Note that here the read alignments with mapping qualities lower than 20 are filtered out. `-Q 20`
`depth2cov` is a binary executable that converts the output of `samtools depth` to `.cov` format.

### 3. Setting Thresholds
The frequencies of coverages can be calculated w/ `cov2counts` which is a binary executable.
````
./cov2counts -i read_alignment.cov -o read_alignment.counts
````
The output file with `.counts` suffix is a 2-cloumn tab-delimited file; the first column shows the coverages and the second column shows the frequencies of
those coverages.
To set the lower- and higher-bound coverage thresholds we calculate the mean and standard deviation of the base-level coverages along the assembly.
We can do this by running the python script `calc_mean_sd.py`.
````
python3 calc_mean_sd.py --countsInput read_alignment.counts --meanOutput cov_mean.txt --sdOutput cov_sd.txt
````
Then we can use the formulas below to calculate the thresholds:

![equation](https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/asset/min_formula.png)

![equation](https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/asset/max_formula.png)

They are actually calculated by these commands:

````
# calculate the min coverage threshold for asset
# using the formula, min( max(5, mean - 3 x sd), 10)
MIN_COVERAGE_ASSET=`awk -v mean=~{coverageMean} -v sd=~{coverageSD} 'BEGIN {min_cov = mean - 3 * sd; if (min_cov < 5) {min_cov=5}; if (min_cov > 10) {min_cov=10}; printf "%d",min_cov}'`

# calculate the min coverage threshold for asset
# using the formula, max(2.5 * mean, mean + 3 x sd)
MAX_COVERAGE_ASSET=`awk -v mean=~{coverageMean} -v sd=~{coverageSD} 'BEGIN {max_cov = mean + 3 * sd; if (max_cov < (2.5 * mean)) {max_cov=2.5 * mean}; printf "%d",max_cov}'`
````
 
 ### 4. Running Asset
 
 Asset is a tool designed for evaluating draft assemblies. It has different modules each of which works with a specfic platform. For more information
 you can visit [Asset github page](https://github.com/dfguan/asset).
 
 Since our platform is PacBio HiFi we used the module `ast_pb`.
 As inputs, it receives an alignment file in `.paf` format and two thresholds that we have already calculated in the previous section.
 
 We converted the bam file we produced in section 1 to `.paf` format with `paftools.js` which is available [here](https://github.com/lh3/minimap2/tree/master/misc). 
 Then we pass the paf file to `ast_pb`:
 ````
 ast_pb -m ${MIN_COVERAGE_ASSET} -M ${MAX_COVERAGE_ASSET} alignment.paf > reliable_blocks.bed
 ````
`reliable_blocks.bed` will contain the blocks whose coverages fall within the given thresholds.

### 5. Post-Processing

We subtracted the whole assembly regions from `realible_blocks.bed` to get the unreliable blocks.
````
# Get the contigs length
samtools faidx asm.fa
 
# Produce a bed file pointing to the whole assembly regions
awk '{print $1"\t0\t"$2}' asm.fa.fai > asm.bed

# Get the unreliable blocks
bedtools subtract -a asm.bed -b reliable_blocks.bed | bedtools merge -d 100 -i - > unreliable_blocks.bed
````
Since the head and tail bases of each contig are usually having low coverages we removed the first and last 1Kb of each contig from the unreliable blocks.
````
# Get contig ends while ignoring contigs <2kb
cat .asm.bed | awk '($3-$2) > 2000 {print $1"\t0\t1000\n"$1"\t"($3-1000)"\t"$3}' > asm.ends.bed

bedtools subtract -a unreliable_blocks -b asm.ends.bed > unreliable_blocks.trim1k.bed
````

### Results Availability

The results are available here: (**These results are for v2 (and v2.1) HPRC_PLUS assemblies)

https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=submissions/8866a160-a3fb-11eb-b821-0a13c5208311--COVERAGE_ANALYSIS/HPRC_PLUS/

As an example here is the list of results for HG002 assemblies:
````
HG002.maternal.hifi.coverage_dist.png
HG002.maternal.hifi.low_high.trim1k.bed
HG002.paternal.hifi.coverage_dist.png
HG002.paternal.hifi.low_high.trim1k.bed
````
It contains one figure and one bed file per haplotype. Each figure shows the coverage distribution and the lower- and higher-bound thresholds are drawn 
with vertical red lines. Each bed is equivalent to `unreliable_blocks.trim1k.bed` similar to what was produced in section 5.

### Codes Availability 

All source codes and scripts are available in `asset/docker/asset/scripts/` or `coverage/docker/coverage/scripts`

### Acknowledgments

Thanks to Arang Rhie for sharing her experience of working with Asset.
