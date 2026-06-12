version 1.0

import "../tasks/compleasm.wdl" as compleasm_wf
import "../tasks/asmgene.wdl" as asmgene_wf
import "../tasks/misjoinCheck.wdl" as misjoin_check_wf
import "../tasks/find_assembly_breakpoints.wdl" as find_breakpoints_wf
import "../tasks/yak_non_trio.wdl" as yak_qv_wf
import "../tasks/yak_no_qv.wdl" as yak_trio_wf
import "../tasks/merqury.wdl" as merqury_wf
import "../tasks/meryl.wdl" as meryl_wf
import "../tasks/extract_reads.wdl" as extract_reads_wf

workflow assembly_qc {

    meta {
        author: "Julian Lucas"
        email: "juklucas@ucsc.edu"
        description: "High-level assembly QC workflow Human Pangenome Reference Consortium assemblies."
    }

    input {
        File hap1_fasta
        File hap2_fasta
        String sample_id

        ## compleasm inputs
        File hap1_compleasm_db
        File hap2_compleasm_db
        String hap1_compleasm_lineage = "primates"
        String hap2_compleasm_lineage = "primates"

        ## asmgene inputs
        File asmgene_genes_fa
        File genes_on_ref_paf

        ## misjoin check inputs
        File misjoin_ref_fa
        File misjoin_ref_centromere_bed

        ## find_assembly_breakpoints inputs
        File breakpoint_reference
        File breakpoint_annotation_bed
        File breakpoint_annotation_sd
        File breakpoint_annotation_censat

        ## kmer inputs
        Array[File]? child_ilmn
        Array[File]? paternal_ilmn
        Array[File]? maternal_ilmn

        File? child_yak
        File? paternal_yak
        File? maternal_yak

        File? read_extraction_reference_fasta

        Boolean run_yak = true
        Boolean run_merqury = true

        ## runtime controls
        Int preemptible = 1
        Int yak_kmer_size = 31
        String yak_genome_size = "3.2g"
        String yak_min_sequence_length = "100k"

        Int merqury_kmer_size = 21
        Int merqury_meryl_thread_count = 32
        Int merqury_meryl_mem_size_gb = 128
        Int merqury_meryl_disk_size_gb = 512
        String read_extraction_docker = "mobinasri/bio_base:v0.2"
        String merqury_docker = "juklucas/hpp_merqury:latest"
        String yak_docker = "juklucas/hpp_yak:latest"
    }

    parameter_meta {
        sample_id: "Sample name used to label outputs and the final assembly_qc CSV."
        hap1_fasta: "Hap1 assembly FASTA."
        hap2_fasta: "Hap2 assembly FASTA."
        child_ilmn: "Optional child Illumina reads, as CRAM/BAM/FASTQ(gz) files. Used to build child Yak and Meryl databases when prebuilt databases are not supplied."
        maternal_ilmn: "Optional maternal Illumina reads, as CRAM/BAM/FASTQ(gz) files. Used to build maternal Yak and Merqury hapmer databases when needed."
        paternal_ilmn: "Optional paternal Illumina reads, as CRAM/BAM/FASTQ(gz) files. Used to build paternal Yak and Merqury hapmer databases when needed."
        read_extraction_reference_fasta: "Reference FASTA for extracting reads from CRAM inputs."
        child_yak: "Optional prebuilt child Yak database for Yak QV. If absent and child_ilmn is provided, the workflow builds one from reads."
        maternal_yak: "Optional prebuilt maternal Yak database for Yak trio phasing metrics. If absent and maternal_ilmn is provided, the workflow builds one from reads."
        paternal_yak: "Optional prebuilt paternal Yak database for Yak trio phasing metrics. If absent and paternal_ilmn is provided, the workflow builds one from reads."
        run_yak: "Run Yak QV and, when both parental kmer sources are available, Yak trio switch/hamming metrics."
        yak_kmer_size: "Optional Yak kmer size used only when building Yak databases from reads. Ignored when prebuilt Yak databases are supplied."
        run_merqury: "Run Merqury QV. If child and both parental reads are provided, also build hapmers for Merqury trio phasing metrics."
        merqury_kmer_size: "Optional Meryl/Merqury kmer size used when building Meryl databases from reads."
        asmgene_genes_fa: "Transcript or gene FASTA used by asmgene."
        genes_on_ref_paf: "Precomputed alignment of asmgene_genes_fa to the reference genome."
        hap1_compleasm_db: "Compleasm lineage database tarball for hap1."
        hap2_compleasm_db: "Compleasm lineage database tarball for hap2."
        hap1_compleasm_lineage: "Compleasm lineage name matching hap1_compleasm_db."
        hap2_compleasm_lineage: "Compleasm lineage name matching hap2_compleasm_db."
        misjoin_ref_fa: "Reference FASTA for minigraph/paftools misjoin detection."
        misjoin_ref_centromere_bed: "Regions (Centromeric) to ignore or treat specially during misjoin detection."
        breakpoint_reference: "Reference FASTA for find_assembly_breakpoints. (usually CHM13v2)"
        breakpoint_annotation_bed: "Genome feature BED used to annotate assembly breakpoints."
        breakpoint_annotation_sd: "Segmental duplication BED used to annotate assembly breakpoints."
        breakpoint_annotation_censat: "CenSat BED used to annotate assembly breakpoints."
    }

    Boolean has_child_kmers = defined(child_yak) || defined(child_ilmn)
    Boolean has_paternal_kmers = defined(paternal_yak) || defined(paternal_ilmn)
    Boolean has_maternal_kmers = defined(maternal_yak) || defined(maternal_ilmn)
    Boolean run_yak_qv = run_yak && has_child_kmers
    Boolean run_yak_trio = run_yak && has_paternal_kmers && has_maternal_kmers
    Boolean run_merqury_qv = run_merqury && defined(child_ilmn)
    Boolean run_merqury_trio = run_merqury_qv && defined(paternal_ilmn) && defined(maternal_ilmn)
    Boolean need_child_read_extract = defined(child_ilmn) && ((run_yak_qv && !defined(child_yak)) || run_merqury_qv)
    Boolean need_paternal_read_extract = defined(paternal_ilmn) && ((run_yak_trio && !defined(paternal_yak)) || run_merqury_trio)
    Boolean need_maternal_read_extract = defined(maternal_ilmn) && ((run_yak_trio && !defined(maternal_yak)) || run_merqury_trio)

    ## Length stats from short_qc
    call runAssemblyStats as hap1_assembly_stats {
        input:
            assembly = hap1_fasta,
            haplotype = "hap1",
            sample_id = sample_id
    }

    call runAssemblyStats as hap2_assembly_stats {
        input:
            assembly = hap2_fasta,
            haplotype = "hap2",
            sample_id = sample_id
    }

    ## Compleasm
    call compleasm_wf.compleasm as compleasm_hap1 {
        input:
            assembly = hap1_fasta,
            lineage_tar = hap1_compleasm_db,
            lineage = hap1_compleasm_lineage
    }

    call compleasm_wf.compleasm as compleasm_hap2 {
        input:
            assembly = hap2_fasta,
            lineage_tar = hap2_compleasm_db,
            lineage = hap2_compleasm_lineage
    }

    ## asmgene
    call asmgene_wf.asmgene as asmgene_hap1 {
        input:
            assemblyFasta = hap1_fasta,
            genesFasta = asmgene_genes_fa,
            genesToReferencePaf = genes_on_ref_paf
    }

    call asmgene_wf.asmgene as asmgene_hap2 {
        input:
            assemblyFasta = hap2_fasta,
            genesFasta = asmgene_genes_fa,
            genesToReferencePaf = genes_on_ref_paf
    }

    ## Misjoin check
    call misjoin_check_wf.misjoinCheck as misjoin_check_hap1 {
        input:
            asm_fasta = hap1_fasta,
            ref_fasta = misjoin_ref_fa,
            ref_centromere_bed = misjoin_ref_centromere_bed,
            name = "~{sample_id}_hap1"
    }

    call misjoin_check_wf.misjoinCheck as misjoin_check_hap2 {
        input:
            asm_fasta = hap2_fasta,
            ref_fasta = misjoin_ref_fa,
            ref_centromere_bed = misjoin_ref_centromere_bed,
            name = "~{sample_id}_hap2"
    }

    ## Breakpoint analysis
    call find_breakpoints_wf.findAssemblyBreakpoints as find_breakpoints_hap1 {
        input:
            assembly = hap1_fasta,
            assembly_name = "~{sample_id}_hap1",
            reference = breakpoint_reference,
            annotationBed = breakpoint_annotation_bed,
            annotationSD = breakpoint_annotation_sd,
            annotationCENSAT = breakpoint_annotation_censat,
            preemptible = preemptible
    }

    call find_breakpoints_wf.findAssemblyBreakpoints as find_breakpoints_hap2 {
        input:
            assembly = hap2_fasta,
            assembly_name = "~{sample_id}_hap2",
            reference = breakpoint_reference,
            annotationBed = breakpoint_annotation_bed,
            annotationSD = breakpoint_annotation_sd,
            annotationCENSAT = breakpoint_annotation_censat,
            preemptible = preemptible
    }

    ## Extract short reads once for Yak and Meryl.
    if (need_child_read_extract) {
        scatter (read_file in select_first([child_ilmn])) {
            call extract_reads_wf.extractReads as child_reads_for_kmers {
                input:
                    readFile = read_file,
                    referenceFasta = read_extraction_reference_fasta,
                    memSizeGB = 4,
                    threadCount = 4,
                    diskSizeGB = 256,
                    dockerImage = read_extraction_docker
            }
        }
    }

    if (need_paternal_read_extract) {
        scatter (read_file in select_first([paternal_ilmn])) {
            call extract_reads_wf.extractReads as paternal_reads_for_kmers {
                input:
                    readFile = read_file,
                    referenceFasta = read_extraction_reference_fasta,
                    memSizeGB = 4,
                    threadCount = 4,
                    diskSizeGB = 256,
                    dockerImage = read_extraction_docker
            }
        }
    }

    if (need_maternal_read_extract) {
        scatter (read_file in select_first([maternal_ilmn])) {
            call extract_reads_wf.extractReads as maternal_reads_for_kmers {
                input:
                    readFile = read_file,
                    referenceFasta = read_extraction_reference_fasta,
                    memSizeGB = 4,
                    threadCount = 4,
                    diskSizeGB = 256,
                    dockerImage = read_extraction_docker
            }
        }
    }

    ## Build missing Yak DBs from reads.
    if (run_yak_qv && !defined(child_yak)) {
        call yakCountFromReads as child_yak_count {
            input:
                readFiles = select_first([child_reads_for_kmers.extractedRead]),
                sampleName = sample_id,
                kmerSize = yak_kmer_size,
                dockerImage = yak_docker
        }
    }

    if (run_yak_trio && !defined(paternal_yak)) {
        call yakCountFromReads as paternal_yak_count {
            input:
                readFiles = select_first([paternal_reads_for_kmers.extractedRead]),
                sampleName = "~{sample_id}_paternal",
                kmerSize = yak_kmer_size,
                dockerImage = yak_docker
        }
    }

    if (run_yak_trio && !defined(maternal_yak)) {
        call yakCountFromReads as maternal_yak_count {
            input:
                readFiles = select_first([maternal_reads_for_kmers.extractedRead]),
                sampleName = "~{sample_id}_maternal",
                kmerSize = yak_kmer_size,
                dockerImage = yak_docker
        }
    }

    ## Yak QV and trio phasing metrics run independently.
    if (run_yak_qv) {
        call yak_qv_wf.yakNonTrioAssemblyStats as yak_qv {
            input:
                assemblyFastaHap1 = hap1_fasta,
                assemblyFastaHap2 = hap2_fasta,
                sampleYak = select_first([child_yak, child_yak_count.outputYak]),
                genomeSize = yak_genome_size,
                minSequenceLength = yak_min_sequence_length,
                dockerImage = yak_docker
        }
    }

    if (run_yak_trio) {
        call yak_trio_wf.yakAssemblyStats as yak_trio {
            input:
                assemblyFastaPat = hap1_fasta,
                assemblyFastaMat = hap2_fasta,
                patYak = select_first([paternal_yak, paternal_yak_count.outputYak]),
                matYak = select_first([maternal_yak, maternal_yak_count.outputYak]),
                dockerImage = yak_docker
        }
    }

    ## Build Meryl DBs for Merqury
    if (run_merqury_qv) {
        call merylCountFromReads as merqury_sample_meryl {
            input:
                readFiles = select_first([child_reads_for_kmers.extractedRead]),
                identifier = "sample",
                kmerSize = merqury_kmer_size,
                memSizeGB = merqury_meryl_mem_size_gb,
                threadCount = merqury_meryl_thread_count,
                diskSizeGB = merqury_meryl_disk_size_gb,
                dockerImage = merqury_docker
        }
    }

    if (run_merqury_trio) {
        call merylCountFromReads as merqury_paternal_meryl {
            input:
                readFiles = select_first([paternal_reads_for_kmers.extractedRead]),
                identifier = "paternal",
                kmerSize = merqury_kmer_size,
                memSizeGB = merqury_meryl_mem_size_gb,
                threadCount = merqury_meryl_thread_count,
                diskSizeGB = merqury_meryl_disk_size_gb,
                dockerImage = merqury_docker
        }

        call merylCountFromReads as merqury_maternal_meryl {
            input:
                readFiles = select_first([maternal_reads_for_kmers.extractedRead]),
                identifier = "maternal",
                kmerSize = merqury_kmer_size,
                memSizeGB = merqury_meryl_mem_size_gb,
                threadCount = merqury_meryl_thread_count,
                diskSizeGB = merqury_meryl_disk_size_gb,
                dockerImage = merqury_docker
        }

        call meryl_wf.merylHapmer as merqury_hapmers {
            input:
                sampleMerylDB = select_first([merqury_sample_meryl.merylDb]),
                paternalMerylDB = merqury_paternal_meryl.merylDb,
                maternalMerylDB = merqury_maternal_meryl.merylDb,
                dockerImage = merqury_docker
        }
    }

    if (run_merqury_qv) {
        call merqury_wf.merqury as merqury {
            input:
                assemblyFasta = hap1_fasta,
                altHapFasta = hap2_fasta,
                kmerTarball = select_first([merqury_sample_meryl.merylDb]),
                patKmerTarball = merqury_hapmers.paternalHapmers,
                matKmerTarball = merqury_hapmers.maternalHapmers,
                dockerImage = merqury_docker
        }
    }

    call collateAssemblyQcCsv {
        input:
            sample_id = sample_id,
            hap1_len_stats = hap1_assembly_stats.lenStats,
            hap2_len_stats = hap2_assembly_stats.lenStats,
            hap1_asmgene_stats = asmgene_hap1.geneStats,
            hap2_asmgene_stats = asmgene_hap2.geneStats,
            hap1_compleasm_summary = compleasm_hap1.summary,
            hap2_compleasm_summary = compleasm_hap2.summary,
            hap1_misjoin_summary = misjoin_check_hap1.misjoinSummary,
            hap2_misjoin_summary = misjoin_check_hap2.misjoinSummary,
            hap1_t2t_contigs = find_breakpoints_hap1.T2Tcontigs,
            hap2_t2t_contigs = find_breakpoints_hap2.T2Tcontigs,
            hap1_t2t_scaffolds = find_breakpoints_hap1.T2Tscaffolds,
            hap2_t2t_scaffolds = find_breakpoints_hap2.T2Tscaffolds,
            yak_qv_summary = yak_qv.outputSummary,
            yak_trio_summary = yak_trio.outputSummary,
            merqury_qv = merqury.QV,
            merqury_results = merqury.outputTarball
    }

    output {
        File assembly_qc_csv = collateAssemblyQcCsv.assemblyQcCsv

        ## assembly stats
        File hap1_len_stats = hap1_assembly_stats.lenStats
        File hap2_len_stats = hap2_assembly_stats.lenStats

        ## compleasm
        File hap1_compleasm_summary = compleasm_hap1.summary
        File hap1_compleasm_full_table = compleasm_hap1.fullTable
        File hap1_compleasm_tar = compleasm_hap1.outputTar
        File hap2_compleasm_summary = compleasm_hap2.summary
        File hap2_compleasm_full_table = compleasm_hap2.fullTable
        File hap2_compleasm_tar = compleasm_hap2.outputTar

        ## asmgene
        File hap1_asmgene_stats = asmgene_hap1.geneStats
        File hap2_asmgene_stats = asmgene_hap2.geneStats

        ## misjoin check
        File hap1_misjoin_summary = misjoin_check_hap1.misjoinSummary
        File hap2_misjoin_summary = misjoin_check_hap2.misjoinSummary

        ## find_assembly_breakpoints
        File hap1_T2Tcontigs = find_breakpoints_hap1.T2Tcontigs
        File hap1_T2Tscaffolds = find_breakpoints_hap1.T2Tscaffolds
        # File hap1_unifiedAssembly = find_breakpoints_hap1.unifiedAssembly
        # File hap1_breakAnnotation_region = find_breakpoints_hap1.breakAnnotation_region
        # File hap1_breakAnnotation_SD = find_breakpoints_hap1.breakAnnotation_SD
        # File hap1_breakAnnotation_CENSAT = find_breakpoints_hap1.breakAnnotation_CENSAT
        # File hap1_assembly_CHM13 = find_breakpoints_hap1.assembly_CHM13
        # File hap1_filteredFlanksBed = find_breakpoints_hap1.filteredFlanksBed
        # File hap1_bed_region = find_breakpoints_hap1.bed_region
        # File hap1_bed_SD = find_breakpoints_hap1.bed_SD
        # File hap1_bed_CENSAT = find_breakpoints_hap1.bed_CENSAT
        # File hap1_assemblyStatistics = find_breakpoints_hap1.assemblyStatistics

        File hap2_T2Tcontigs = find_breakpoints_hap2.T2Tcontigs
        File hap2_T2Tscaffolds = find_breakpoints_hap2.T2Tscaffolds
        # File hap2_unifiedAssembly = find_breakpoints_hap2.unifiedAssembly
        # File hap2_breakAnnotation_region = find_breakpoints_hap2.breakAnnotation_region
        # File hap2_breakAnnotation_SD = find_breakpoints_hap2.breakAnnotation_SD
        # File hap2_breakAnnotation_CENSAT = find_breakpoints_hap2.breakAnnotation_CENSAT
        # File hap2_assembly_CHM13 = find_breakpoints_hap2.assembly_CHM13
        # File hap2_filteredFlanksBed = find_breakpoints_hap2.filteredFlanksBed
        # File hap2_bed_region = find_breakpoints_hap2.bed_region
        # File hap2_bed_SD = find_breakpoints_hap2.bed_SD
        # File hap2_bed_CENSAT = find_breakpoints_hap2.bed_CENSAT
        # File hap2_assemblyStatistics = find_breakpoints_hap2.assemblyStatistics

        ## Yak outputs
        File? created_child_yak = child_yak_count.outputYak
        File? created_paternal_yak = paternal_yak_count.outputYak
        File? created_maternal_yak = maternal_yak_count.outputYak
        File? yak_qv_summary = yak_qv.outputSummary
        File? yak_qv_tar = yak_qv.outputTarball
        File? yak_trio_summary = yak_trio.outputSummary
        File? yak_trio_tar = yak_trio.outputTarball

        ## Merqury / Meryl outputs
        File? sample_meryl_db = merqury_sample_meryl.merylDb
        File? paternal_meryl_db = merqury_paternal_meryl.merylDb
        File? maternal_meryl_db = merqury_maternal_meryl.merylDb
        File? merqury_maternal_hapmers = merqury_hapmers.maternalHapmers
        File? merqury_paternal_hapmers = merqury_hapmers.paternalHapmers
        File? merqury_hapmer_images = merqury_hapmers.hapmerImages
        File? merqury_qv = merqury.QV
        File? merqury_tar = merqury.outputTarball
    }
}

task runAssemblyStats {
    input {
        File assembly
        String haplotype
        String sample_id

        Int memSizeGB = 16
        Int threadCount = 4
        Int diskSizeGB = 64
        String dockerImage = "humanpangenomics/hpp_qc_stats@sha256:6a64ac0be88ce9ca760eb7713922f65e66f9d09076b9d17c2c416e5d558bf1d0"
    }

    String output_name = "~{sample_id}.~{haplotype}.len.stats.txt"

    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        k8 ${CAL_N50_PATH} ~{assembly} > ~{output_name}
    >>>

    output {
        File lenStats = output_name
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}

task yakCountFromReads {
    input {
        Array[File] readFiles
        String sampleName
        Int kmerSize = 31
        Int bloomSize = 37

        Int memSizeGB = 128
        Int threadCount = 16
        Int diskSizeGB = 256
        String dockerImage = "juklucas/hpp_yak:latest"
    }

    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        yak count -k~{kmerSize} -t~{threadCount} -b~{bloomSize} -o ~{sampleName}.yak <(cat ~{sep=" " readFiles}) <(cat ~{sep=" " readFiles})
    >>>

    output {
        File outputYak = "~{sampleName}.yak"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}

task merylCountFromReads {
    input {
        Array[File] readFiles
        String identifier
        Int kmerSize = 21
        Boolean compress = false

        Int memSizeGB = 128
        Int threadCount = 32
        Int diskSizeGB = 512
        String dockerImage = "juklucas/hpp_merqury:latest"
    }

    String compress_arg = if compress then "compress" else ""

    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace
        export OMP_NUM_THREADS=~{threadCount}

        meryl ~{compress_arg} k=~{kmerSize} threads=~{threadCount} memory=$((~{memSizeGB}-10)) count output ~{identifier}.meryl ~{sep=" " readFiles}

        tar cvf ~{identifier}.meryl.tar ~{identifier}.meryl
        rm -rf ~{identifier}.meryl
    >>>

    output {
        File merylDb = identifier + ".meryl.tar"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}

task collateAssemblyQcCsv {
    input {
        String sample_id

        File hap1_len_stats
        File hap2_len_stats
        File hap1_asmgene_stats
        File hap2_asmgene_stats
        File hap1_compleasm_summary
        File hap2_compleasm_summary
        File hap1_misjoin_summary
        File hap2_misjoin_summary
        File hap1_t2t_contigs
        File hap2_t2t_contigs
        File hap1_t2t_scaffolds
        File hap2_t2t_scaffolds

        File? yak_qv_summary
        File? yak_trio_summary
        File? merqury_qv
        File? merqury_results

        Int memSizeGB = 8
        Int diskSizeGB = 64
        String dockerImage = "python:3.11-slim"
    }

    String output_name = "~{sample_id}.assembly_qc.csv"

    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        python3 <<'PY'
        import csv
        import os
        import re
        import tarfile

        NA = "NA"

        def clean(v):
            if v is None or v == "":
                return NA
            return str(v)

        def first_number(tokens):
            for token in tokens:
                token = token.strip().strip(",;%")
                if re.fullmatch(r"[-+]?[0-9]*\.?[0-9]+", token):
                    return token
            return None

        def count_data_lines(path):
            try:
                with open(path) as fh:
                    return str(sum(1 for line in fh if line.strip() and not line.lstrip().startswith("#") and not line.lstrip().lower().startswith("name")))
            except Exception:
                return NA

        def parse_len_stats(path):
            out = {
                "assembly_size": NA,
                "number_of_contigs": NA,
                "n50": NA,
                "aun": NA,
            }

            def looks_numeric(value):
                return re.fullmatch(r"[-+]?[0-9][0-9,]*(\.[0-9]+)?", str(value).strip()) is not None

            try:
                lines = [line.rstrip("\n") for line in open(path)]
                if len(lines) > 19:
                    indexed = {
                        "assembly_size": lines[6].split("\t")[1],
                        "number_of_contigs": lines[7].split("\t")[1],
                        "n50": lines[13].split("\t")[2],
                        "aun": lines[19].split("\t")[1],
                    }
                    if all(looks_numeric(value) for value in indexed.values()):
                        out.update(indexed)
                        return out
            except Exception:
                pass

            try:
                for line in open(path):
                    fields = re.split(r"\s+", line.strip())
                    low = line.lower()
                    if "total" in low and "length" in low:
                        out["assembly_size"] = clean(first_number(fields[1:]))
                    elif "contig" in low and ("number" in low or "count" in low):
                        out["number_of_contigs"] = clean(first_number(fields[1:]))
                    elif "n50" in low:
                        out["n50"] = clean(first_number(fields[1:]))
                    elif "aun" in low:
                        out["aun"] = clean(first_number(fields[1:]))
            except Exception:
                pass
            return out

        def parse_asmgene(path):
            labels = [
                "full_sgl",
                "full_dup",
                "frag_sgl",
                "frag_dup",
                "part_sgl",
                "part_dup",
                "missing",
            ]
            out = {f"{label}_asmgene": NA for label in labels}
            try:
                for line in open(path):
                    fields = re.split(r"\s+", line.strip())
                    normalized_fields = [field.strip().rstrip(":") for field in fields]
                    for label in labels:
                        match = re.search(rf"\b{re.escape(label)}\b[:\s]+([-+]?[0-9]*\.?[0-9]+)", line)
                        if match:
                            out[f"{label}_asmgene"] = clean(match.group(1))
                        elif label in normalized_fields:
                            i = normalized_fields.index(label)
                            out[f"{label}_asmgene"] = clean(first_number(fields[i + 1:]))
            except Exception:
                pass
            return out

        def parse_compleasm(path):
            out = {
                "sgl_compleasm": NA,
                "dup_compleasm": NA,
                "fragmented_compleasm": NA,
                "missing_compleasm": NA,
                "total_compleasm": NA,
            }
            code_to_metric = {
                "S": "sgl_compleasm",
                "D": "dup_compleasm",
                "F": "fragmented_compleasm",
                "M": "missing_compleasm",
                "N": "total_compleasm",
            }
            try:
                for line in open(path):
                    m = re.match(r"\s*([SDFMN]):", line)
                    if not m:
                        continue
                    nums = re.findall(r"[-+]?[0-9]*\.?[0-9]+", line)
                    if nums:
                        out[code_to_metric[m.group(1)]] = nums[-1]
            except Exception:
                pass
            return out

        def parse_yak_qv(path):
            out = {"hap1": NA, "hap2": NA}
            if not path or not os.path.exists(path):
                return out
            current = None
            last_data = None
            try:
                for line in open(path):
                    stripped = line.strip()
                    low = stripped.lower()
                    if low.startswith("#"):
                        if current and last_data:
                            out[current] = clean(last_data.split()[-1])
                        current = None
                        last_data = None
                        if "hap1" in low or "pat qv" in low:
                            current = "hap1"
                        elif "hap2" in low or "mat qv" in low:
                            current = "hap2"
                    elif current and stripped:
                        last_data = stripped
                if current and last_data:
                    out[current] = clean(last_data.split()[-1])
            except Exception:
                pass
            return out

        def parse_yak_trio(path):
            out = {
                "switch_yak": {"hap1": NA, "hap2": NA},
                "hamming_yak": {"hap1": NA, "hap2": NA},
            }
            if not path or not os.path.exists(path):
                return out
            current = None
            data = []

            def flush():
                if not current or not data:
                    return
                hap = current
                if len(data) >= 1:
                    out["switch_yak"][hap] = clean(data[0].split()[-1])
                if len(data) >= 2:
                    out["hamming_yak"][hap] = clean(data[1].split()[-1])

            try:
                for line in open(path):
                    stripped = line.strip()
                    low = stripped.lower()
                    if low.startswith("#"):
                        flush()
                        data = []
                        current = None
                        if "pat" in low or "hap1" in low:
                            current = "hap1"
                        elif "mat" in low or "hap2" in low:
                            current = "hap2"
                    elif current and stripped:
                        data.append(stripped)
                flush()
            except Exception:
                pass
            return out

        def parse_merqury_qv(path):
            out = {"hap1": NA, "hap2": NA}
            if not path or not os.path.exists(path):
                return out
            try:
                for line in open(path):
                    fields = line.strip().split()
                    low = line.lower()
                    if not fields:
                        continue
                    value = fields[3] if len(fields) > 3 else fields[-1]
                    if "althap" in low:
                        out["hap2"] = value
                    elif re.search(r"(^|\s)asm(\s|$)", low):
                        out["hap1"] = value
            except Exception:
                pass
            return out

        def parse_merqury_phase(path):
            out = {
                "switch_merqury": {"hap1": NA, "hap2": NA},
            }
            if not path or not os.path.exists(path):
                return out
            extract_dir = "merqury_extract"
            os.makedirs(extract_dir, exist_ok=True)
            try:
                with tarfile.open(path) as tar:
                    extract_root = os.path.abspath(extract_dir)
                    for member in tar.getmembers():
                        target = os.path.abspath(os.path.join(extract_root, member.name))
                        if target == extract_root or target.startswith(extract_root + os.sep):
                            tar.extract(member, extract_root)
            except Exception:
                return out

            for root, _, files in os.walk(extract_dir):
                for name in files:
                    fn = os.path.join(root, name)
                    try:
                        with open(fn, errors="ignore") as fh:
                            for line in fh:
                                low = line.lower()
                                if "switch" not in low:
                                    continue
                                nums = re.findall(r"[-+]?[0-9]*\.?[0-9]+", line)
                                if not nums:
                                    continue
                                value = nums[-1]
                                if "althap" in low or "hap2" in low or ".mat" in low:
                                    out["switch_merqury"]["hap2"] = value
                                elif "asm" in low or "hap1" in low or ".pat" in low:
                                    out["switch_merqury"]["hap1"] = value
                    except Exception:
                        continue
            return out

        def parse_misjoins(path):
            try:
                count = 0
                for line in open(path):
                    stripped = line.strip()
                    if not stripped or stripped.startswith("#"):
                        continue
                    fields = stripped.split()
                    # Only uppercase J records are inter-chromosomal misjoins; other records can reflect mismapping elsewhere.
                    if fields and fields[0] == "J":
                        count += 1
                return str(count)
            except Exception:
                return NA

        def write_row(writer, metric, hap1, hap2):
            writer.writerow([metric, clean(hap1), clean(hap2)])

        len_hap1 = parse_len_stats("~{hap1_len_stats}")
        len_hap2 = parse_len_stats("~{hap2_len_stats}")
        asmgene_hap1 = parse_asmgene("~{hap1_asmgene_stats}")
        asmgene_hap2 = parse_asmgene("~{hap2_asmgene_stats}")
        compleasm_hap1 = parse_compleasm("~{hap1_compleasm_summary}")
        compleasm_hap2 = parse_compleasm("~{hap2_compleasm_summary}")
        yak_qv = parse_yak_qv("~{yak_qv_summary}")
        yak_trio = parse_yak_trio("~{yak_trio_summary}")
        merqury_qv = parse_merqury_qv("~{merqury_qv}")
        merqury_phase = parse_merqury_phase("~{merqury_results}")

        t2t_sequences_hap1 = count_data_lines("~{hap1_t2t_scaffolds}")
        t2t_sequences_hap2 = count_data_lines("~{hap2_t2t_scaffolds}")
        t2t_contigs_hap1 = count_data_lines("~{hap1_t2t_contigs}")
        t2t_contigs_hap2 = count_data_lines("~{hap2_t2t_contigs}")

        def subtract(a, b):
            try:
                return str(max(int(a) - int(b), 0))
            except Exception:
                return NA

        with open("~{output_name}", "w", newline="") as out:
            writer = csv.writer(out)
            writer.writerow(["metric", "hap1", "hap2"])

            write_row(writer, "qv_yak", yak_qv["hap1"], yak_qv["hap2"])
            write_row(writer, "switch_yak", yak_trio["switch_yak"]["hap1"], yak_trio["switch_yak"]["hap2"])
            write_row(writer, "hamming_yak", yak_trio["hamming_yak"]["hap1"], yak_trio["hamming_yak"]["hap2"])
            write_row(writer, "qv_merqury", merqury_qv["hap1"], merqury_qv["hap2"])
            write_row(writer, "switch_merqury", merqury_phase["switch_merqury"]["hap1"], merqury_phase["switch_merqury"]["hap2"])

            for metric in [
                "full_sgl_asmgene",
                "full_dup_asmgene",
            ]:
                write_row(writer, metric, asmgene_hap1[metric], asmgene_hap2[metric])

            for metric in [
                "sgl_compleasm",
                "dup_compleasm",
                "fragmented_compleasm",
                "missing_compleasm",
                "total_compleasm",
            ]:
                write_row(writer, metric, compleasm_hap1[metric], compleasm_hap2[metric])

            write_row(writer, "assembly_size", len_hap1["assembly_size"], len_hap2["assembly_size"])
            write_row(writer, "number_of_contigs", len_hap1["number_of_contigs"], len_hap2["number_of_contigs"])
            write_row(writer, "n50", len_hap1["n50"], len_hap2["n50"])
            write_row(writer, "aun", len_hap1["aun"], len_hap2["aun"])
            write_row(writer, "t2t_sequences", t2t_sequences_hap1, t2t_sequences_hap2)
            write_row(writer, "t2t_contigs", t2t_contigs_hap1, t2t_contigs_hap2)
            write_row(writer, "t2t_scaffolds", subtract(t2t_sequences_hap1, t2t_contigs_hap1), subtract(t2t_sequences_hap2, t2t_contigs_hap2))
            write_row(writer, "misjoins", parse_misjoins("~{hap1_misjoin_summary}"), parse_misjoins("~{hap2_misjoin_summary}"))
        PY
    >>>

    output {
        File assemblyQcCsv = output_name
    }

    runtime {
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
