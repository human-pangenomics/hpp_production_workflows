version 1.0

import "../tasks/hifiasm_cutadapt_multistep.wdl" as hifiasm_multistep_wf
import "../tasks/gfatools.wdl" as gfatools_wf

workflow trioHifiasmAssembly {
    input {
        String childID
        Array[File] childReadsHiFi
        Array[File] childReadsONT=[]
        Int? minOntReadLength
        Int? homCov
        File paternalYak
        File maternalYak
        Boolean filterAdapters=true
        String? hifiasmExtraOptions
        File? inputBinFilesTarGz
        File? referenceFasta
        # runtime configurations for running hifiasm
        Int threadCount=64
        Int preemptible=2
    }

    

    ### Trio Hifiasm ###
    call hifiasm_multistep_wf.runTrioHifiasm as trioHifiasm{
        input:
            paternalYak = paternalYak,
            maternalYak = maternalYak,
            childReadsHiFi = childReadsHiFi,
            childReadsONT = childReadsONT,
            minOntReadLength = minOntReadLength,
            homCov = homCov,
            childID = childID,
            filterAdapters = filterAdapters,
            hifiasmExtraOptions = hifiasmExtraOptions,
            inputBinFilesTarGz = inputBinFilesTarGz,
            referenceFasta = referenceFasta,
            threadCount = threadCount,
            preemptible = preemptible
    }

    ### Convert GFA to FASTA ###
    call gfatools_wf.phasedGFAs2Fasta as gfa2fasta{
        input:
            paternalGfa = trioHifiasm.outputPaternalGfa,
            maternalGfa = trioHifiasm.outputMaternalGfa,
            childID = childID
    }

    output {
        File paternalFastaGz = gfa2fasta.outputPaternalFastaGz
        File maternalFastaGz = gfa2fasta.outputMaternalFastaGz
        File paternalContigGfaTarGz = trioHifiasm.outputPaternalContigGfa 
        File maternalContigGfaTarGz = trioHifiasm.outputMaternalContigGfa 
        File rawUnitigGfaTarGz = trioHifiasm.outputRawUnitigGfa
        File binFilesTarGz = trioHifiasm.outputBinFiles
    }
    parameter_meta {
        childID: "Sample ID of the child sample whose reads are going to be assembled"
        childReadsHiFi: "An array of files (or a single file) that contain the HiFi reads of the child sample ( Acceptable formats are fastq (or fq), fastq.gz (or fq.gz), bam and cram)"
        inputBinFilesTarGz: "(optional) The hifiasm produces some bin files which can be saved and used for re-running the assembly process. By having these bin files hifiasm can skip the time-consuming process of finding overlaps (Acceptable format is .tar.gz)"
        referenceFasta: "(optional) If any of the read files (can be either for child, father or mother) are having .cram format, the reference genome should be provided in .fasta format"
        threadCount: "The number of cores for running hifiasm"
        memSize: "The memory size (GB) for running hifiasm"
        preemptible: "The number of tries for using a preemptible node for running hifiasm. Note that if your child data has a coverage of more than 40X, hifiasm (without any given bin files) may take longer than 24 hours. So using a preemptible node is useless beacuse it gets interrupted after 24 hours"
    }
    meta {
        author: "Mobin Asri"
        email: "masri@ucsc.edu"
    }
}

