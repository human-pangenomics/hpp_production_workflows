version 1.0

workflow fast5GuppyGPU {
    input {
        File fast5_folder_path
        String sample_name
        String guppy_version
        Int desired_size_GB
        
    }

    parameter_meta {
        fast5_list: "List of gs path to fast5 files"
        sample_name: "Name of sample, used for output file names"
        guppy_version: "Guppy version, used for output file names"
        desired_size_GB: "Choose size to split input tar file by. With a 300GB fast5_tar_file and 30GB desired_size_GB, the fast5_tar_file will be split in 10 pieces."
    }

    call pathToList {
        input:
            folder_path = fast5_folder_path
    }

    Array[File] fast5_paths = read_lines(pathToList.path_list)


    call splitFast5s {
        input:
            files_to_split = fast5_paths,
            desired_size_GB = desired_size_GB
    }

    # call guppyGPU on each of the smaller "split" tar files
    scatter (split_fast5 in splitFast5s.split_fast5s) {
        call guppyGPU {
            input:
                fast5_tar_file = split_fast5
        }
    }

    call concatenateBam as passBam {
        input:
            files = flatten(guppyGPU.pass_bam),
            sample_name = sample_name,
            guppy_version = guppy_version,
            pass_fail = "pass"
    }

    call concatenateFastq as passFastq {
        input:
            files = flatten(guppyGPU.pass_fastq),
            sample_name = sample_name,
            guppy_version = guppy_version,
            pass_fail = "pass"
    }

    call concatenateBam as failBam {
        input:
            files = flatten(guppyGPU.fail_bam),
            sample_name = sample_name,
            guppy_version = guppy_version,
            pass_fail = "fail"
    }

    call concatenateFastq as failFastq {
        input:
            files = flatten(guppyGPU.fail_fastq),
            sample_name = sample_name,
            guppy_version = guppy_version,
            pass_fail = "fail"
    }


    call concatenateSummary {
        input:
            files = guppyGPU.summary,
            sample_name = sample_name,
            guppy_version = guppy_version
    }


    output {
        File bams_pass = passBam.concatenatedBam
        File fastqs_pass = passFastq.concatenatedFastq
        File bams_fail = failBam.concatenatedBam
        File fastqs_fail = failFastq.concatenatedFastq
        File summaries = concatenateSummary.concatenatedSummary
    }
        
    meta {
        author: "Jimin Park"
        email: "jpark621@ucsc.edu"
        description: "Calls guppy_basecaller with GPUs. Takes in fast5 files and outputs unaligned bam with methylation calls, fastq and sequencing summary text file."
    }

}

task pathToList {
    input {
        String folder_path
        String dockerImage = "google/cloud-sdk:latest"
        Int preempts = 3
        Int memSizeGB = 8
        Int extraDisk = 5
        Int threadCount = 2
        Int diskSizeGB = 5
    }

    command <<<
        gsutil ls ~{folder_path} > "fast5.list"
    >>>

    output {
        File path_list = "fast5.list"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible : preempts
    }
}


task splitFast5s {
    input {
        Array[File] files_to_split
        Int desired_size_GB

        String dockerImage = "jiminpark/guppy-wdl:latest" 

        Int preempts = 3
        Int memSizeGB = 8
        Int extraDisk = 5
        Int threadCount = 2
    }

    Int file_size = ceil(size(files_to_split, "GB"))
    Int diskSizeGB = 3 * file_size + extraDisk

    command <<<
        # move files into folder until exceeds desired_size_GB
        # then tar contents of folder
        OUTPUT_IDX=0
        OUTPUT_DIR=fast5_tar_$OUTPUT_IDX
        mkdir $OUTPUT_DIR
        for FILE in ~{sep=' ' files_to_split}
        do
            size=$(du -s -BG $OUTPUT_DIR | sed 's/G.*//')
            if (( $size > ~{desired_size_GB} ))
            then
                tar -cvf fast5_tarball_$OUTPUT_IDX.tar $OUTPUT_DIR/*
                rm -r $OUTPUT_DIR
                OUTPUT_IDX=$(($OUTPUT_IDX + 1))
                OUTPUT_DIR=fast5_tar_$OUTPUT_IDX
                mkdir $OUTPUT_DIR
            fi
            mv $FILE $OUTPUT_DIR
        done
        tar -cvf fast5_tarball_$OUTPUT_IDX.tar $OUTPUT_DIR/*
        rm -r $OUTPUT_DIR

    >>>

    output {
        Array[File] split_fast5s = glob("*tar")
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible : preempts
    }
}


task guppyGPU {
    
    input {
    
        File fast5_tar_file

        String CONFIG_FILE = "dna_r9.4.1_450bps_modbases_5mc_cg_sup_prom.cfg"
        Int READ_BATCH_SIZE = 250000
        Int q = 250000

        String dockerImage = "jiminpark/guppy-wdl:latest" 

        String? additionalArgs

        Int preempts = 3
        Int memSizeGB = 64
        Int threadCount = 12
        Int extraDisk = 5
        Int gpuCount = 1
        Int maxRetries = 4
        String gpuType = "nvidia-tesla-v100"
        String nvidiaDriverVersion = "418.87.00"
        String zones = "us-west1-b"
    }

    # calculate needed disk size
    Int file_size = ceil(size(fast5_tar_file, "GB"))
    Int diskSizeGB = 3 * file_size + extraDisk


    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        ## Extract tar file to 
        mkdir input
        
        # place all extracted files into directory input
        tar xvf "~{fast5_tar_file}" --directory input

        mkdir output

        # check if length of "additionalArgs" is zero

        if [[ "~{additionalArgs}" == "" ]]
        then
            ADDITIONAL_ARGS=""
        else
            ADDITIONAL_ARGS="~{additionalArgs}"
        fi


        guppy_basecaller \
            -i input/ \
            -s output/ \
            -c /opt/ont/guppy/data/"~{CONFIG_FILE}" \
            --bam_out \
            -x cuda:all:100% \
            -r \
            --read_batch_size "~{READ_BATCH_SIZE}" \
            -q "~{q}" \
            ${ADDITIONAL_ARGS}

    >>>

    output {
        Array[File] pass_bam = glob("output/*.bam")
        Array[File] pass_fastq = glob("output/*.fastq")

        Array[File] fail_bam = glob("output/fail/*bam")
        Array[File] fail_fastq = glob("output/fail/*fastq")
        File summary = glob("output/sequencing_summary.txt")[0]
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        gpuCount: gpuCount
        gpuType: gpuType
        maxRetries : maxRetries
        nvidiaDriverVersion: nvidiaDriverVersion
        preemptible : preempts
        docker: dockerImage
        zones: zones
    }
}


task concatenateBam {
    input {
        Array[File] files
        
        String sample_name
        String guppy_version
        String pass_fail

        String dockerImage = "tpesout/megalodon:latest"

        Int preempts = 3
        Int memSizeGB = 8
        Int threadCount = 3
        Int diskSizeGB = 500
    }
    
    command {
        samtools merge -o "${sample_name}_${guppy_version}_${pass_fail}.bam" ${sep=" " files}
    }

    output {
        File concatenatedBam = "${sample_name}_${guppy_version}_${pass_fail}.bam"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible : preempts
    }
}


task concatenateFastq {
    input {
        Array[File] files
        String sample_name
        String guppy_version
        String pass_fail

        String dockerImage = "tpesout/megalodon:latest"

        # runtime
        Int preempts = 3
        Int memSizeGB = 8
        Int threadCount = 3
        Int diskSizeGB = 500
    }
    
    command {
        cat ${sep=" " files} | gzip -c > "${sample_name}_${guppy_version}_${pass_fail}.fastq.gz"
    }

    output {
        File concatenatedFastq = "${sample_name}_${guppy_version}_${pass_fail}.fastq.gz"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible : preempts
    }
}


task concatenateSummary {
    input {
        Array[File] files
        String sample_name
        String guppy_version

        String dockerImage = "tpesout/megalodon:latest"

        # runtime
        Int preempts = 3
        Int memSizeGB = 8
        Int threadCount = 3
        Int diskSizeGB = 50
    }
    
    command {
        cat ${sep=" " files} > "tmp.txt"
        # remove duplicate headers
        awk 'NR==1 || !/^filename/' "tmp.txt" > "${sample_name}_${guppy_version}_sequencing_summary.txt"
    }

    output {
        File concatenatedSummary = "${sample_name}_${guppy_version}_sequencing_summary.txt"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible : preempts
    }
}