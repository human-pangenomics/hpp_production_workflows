## CopEval 

### Overview
**CopEval** is a read-based pipeline that can find different types of mis-assemblies in a draft diploid assembly. One core component of this pipeline is another pipeline named **Partitioner**. Partitioner is a pipeline that recieves the read alignments to the draft assembly and partition the assembly into 4 main components; erroneous, (falsely) duplicated, haploid and collapsed. It is explained [here](https://github.com/human-pangenomics/hpp_production_workflows/edit/asset/coverage/docs/coverage/README.md) in detail.

CopEval has 7 steps:
- Aligning reads to the diploid assembly
- Phase the ambiguous alignments using [the phasing pipeline](https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/docs/phasing/README.md)
- Run Partitioner on the alignments
- Call variants 
- Remove the alignments with alternative alleles
- Run Paritioner on the alignments with no alternative allele
- Combine the Partitioner outputs
### 1. Read Alignment
The ONT and HiFi reads are aligned to each diploid assembly (~ 6Gbases long) with winnowmap v2.03. Since we are aligning the reads to the diploid assembly the expected base-level coverage should be half of the sequencing coverage.
Here are the main commands used for producing the alignments (taken from the [winnowmap docs](https://github.com/marbl/Winnowmap)):
```` 
  # making the k-mer table with meryl
  meryl count k=15 output merylDB asm.fa
  meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt
  
  # alignment with winnowmap
  winnowmap -W repetitive_k15.txt -ax [map-ont | map-pb] -Y -L --eqx --cs -I8g <(cat pat_asm.fa mat_asm.fa) reads.fq.gz | \
    samtools view -hb > read_alignment.bam
````

3. Some regions are falsely flagged as collapsed. The reason is that the equivalent region in the other haplotype is not assembled correctly so the reads from two haplotypes are aligned to only one of them. This flagging can be useful since it points to a region whose counterpart in the other haplotype is not assembled correctly or not assembled at all. 

    One approach is to correct the coverage in the correctly assembled region. The coverage can be corrected by detecting the marker snps and removing the reads from the wrong halpotype or segment. This is going to be incorporated in the next releases of this analysis. Here is an example of a region with ~40X coverage but after detecting the marker snps (by variant calling) and removing the wrong alignments the coverage has decreased to ~17X which is much closer to the expected coverage (~20X).

   <img src="https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/images/coverage_correction.png" width="700" height="275">


### Data, Source Code and Workflows Availability

The haplotype-resolved assemblies of the HPRC-Y1 samples and their corresponding data sets are available in

https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=working/HPRC/

For more details read this github page:

https://github.com/human-pangenomics/HPP_Year1_Assemblies

We have used the Genbank version of the HPRC-Y1 assemblies.

The Python scripts, C source codes and the binary files are available in

https://github.com/human-pangenomics/hpp_production_workflows/tree/asset/coverage/docker/coverage/scripts

The wdl files that have been used for this analysis are available in

https://github.com/human-pangenomics/hpp_production_workflows/tree/asset/coverage/wdl/tasks

### Results Availability

The results are available in 

https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=submissions/e9ad8022-1b30-11ec-ab04-0a13c5208311--COVERAGE_ANALYSIS_Y1_GENBANK/

The most updated results are under the `V1` directory. It includes the results of the HiFi-based coverage analysis.
The directory structure for the results of each sample is as below:

````
HG002
├── HG002.diploid.f1_assembly_v2_genbank.hifi.winnowmap_v2.03.cov_dist.pdf
├── HG002.diploid.f1_assembly_v2_genbank.hifi.winnowmap_v2.03.table
├── HG002.diploid.f1_assembly_v2_genbank.hifi.winnowmap_v2.03.tables.tar.gz
├── combined
├── contig_based
├── filtered
├── dup_corrected
└── whole_genome_based
````

The pdf file contains the plots of the coverage distributions for the whole assembly and the contigs longer than 5Mb which has been used for the contig-based analysis. Different model components are also drawn with different colors.

The table files contain the models fit to the distributions. These models are used for finding the thresholds and categorizing the regions into one of the four components; error, duplicated, haploid and collapsed.

The four subdirectories are having the categorized blocks. They are the results of different parts of the pipeline and named based on this.
For example `combined` directory contains the combination of `contig_based` and `whole_genome_based` blocks (section 6) and `dup_corrected` contains the corrections explained in section 7.
Each directory has four bed files and they are named based on the component they are pointing to.

````
filtered/
├── HG002.diploid.f1_assembly_v2_genbank.hifi.winnowmap_v2.03.filtered.collapsed.bed
├── HG002.diploid.f1_assembly_v2_genbank.hifi.winnowmap_v2.03.filtered.duplicated.bed
├── HG002.diploid.f1_assembly_v2_genbank.hifi.winnowmap_v2.03.filtered.error.bed
└── HG002.diploid.f1_assembly_v2_genbank.hifi.winnowmap_v2.03.filtered.haploid.bed
````

For example the bed file that ends with `.haploid.bed` is pointing to the haploid blocks.

It is recommended to use the bed files in the `filtered` subdirectory.

### Acknowledgements

Thanks to Jordan M Eizenga for sharing the source code of EM that fits the mixture model and his useful support for doing this analysis.


