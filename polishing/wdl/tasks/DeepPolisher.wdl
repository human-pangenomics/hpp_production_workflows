version 1.0

# This is a task level wdl workflow to run DeepPolisher

workflow runDeepPolisher {
    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "Run DeepPolisher on phased reads aligned to one haplotype"
    }
    input {
        File Bam
        File Bai
        File Fasta
        File ModelFilesTarGZ
        String sampleName
        String dockerImage="google/deepconsensus:polisher_v0.0.8_12122023"
        Boolean useOptimalGQFilter=true
        String customGQFilter=""
    }
    call DeepPolisher {
        input:
            Bam=Bam,
            Bai=Bai,
            Fasta=Fasta,
            ModelFilesTarGZ=ModelFilesTarGZ,
            sampleName=sampleName,
            dockerImage=dockerImage
    }
    call DPPostProcess {
        input:
            VCFsTarGz=DeepPolisher.VCFsTarGz,
            useOptimalGQFilter=useOptimalGQFilter,
            customGQFilter=customGQFilter
    }
    output {
        File PolisherVcf = DPPostProcess.vcfFile
        File PolisherVcfTbi = DPPostProcess.vcfFileTbi
        File PolisherVcfNoFilter=DPPostProcess.noFiltersPolisherVcf
    }
}

task DeepPolisher{
    input {
        File Bam
        File Bai
        File Fasta
        File ModelFilesTarGZ
        String sampleName

        String dockerImage
        Int memSizeGB = 128
        Int threadCount = 32
        Int diskSizeGB = 256
    }
    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        mkdir images
        mkdir vcf

        # softlink bam and index so they are in same directory
        BAMFILE=$(basename ~{Bam})
        BAIFILE=$(basename ~{Bai})

        ln -s ~{Bam} ./${BAMFILE}
        ln -s ~{Bai} ./${BAIFILE}

        # untar model files
        # they need to be tar'd in a folder called "checkpoint" and the name of the tar file needs to match
        # the prefix for the .autofdo_profile, .index and .data-00000-of-00001 files example:
        # tar -zcvf checkpoint-657.tar.gz checkpoint/

        CHECKPOINT=( $(basename ~{ModelFilesTarGZ} | sed 's/.gz$//' | sed 's/.tar$//') )

        echo $CHECKPOINT

        echo checkpoint/${CHECKPOINT}

        tar xvf ~{ModelFilesTarGZ} --no-same-owner

        # Make images
        time polisher make_images --bam ${BAMFILE} --fasta ~{Fasta} --output images/images --cpus ~{threadCount}

        # Inference on images to generate VCFs
        time polisher inference --input_dir images --out_dir vcf/ --checkpoint checkpoint/${CHECKPOINT} --reference_file ~{Fasta} --sample_name ~{sampleName} --cpus ~{threadCount}

        tar -zcvf vcf.tar.gz vcf/

    >>>
    output {
        File VCFsTarGz="vcf.tar.gz"
    }
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}

task DPPostProcess{
    input{
        File VCFsTarGz

        String customGQFilter=""
        Boolean useOptimalGQFilter=true
        String dockerImage = "miramastoras/polishing:latest"
        Int memSizeGB = 128
        Int threadCount = 8
        Int diskSizeGB = 128
    }
    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        # untar vcf files
        tar xvf ~{VCFsTarGz} --no-same-owner

        # get vcf names
        VCF_NAMES=$(ls vcf/*vcf.gz | tr '\n' ' ')

        # combine vcfs
        vcf-concat $VCF_NAMES > polisher_output.unsorted.vcf

        # Sort the calls
        bcftools view polisher_output.unsorted.vcf --no-header | vcf-sort > polisher_output.sorted.calls_only

        # Get the header
        bcftools view polisher_output.unsorted.vcf --header > polisher_output.header

        # remove stars for now until bug is fixed
        grep -v "*" polisher_output.sorted.calls_only > tmp; mv tmp polisher_output.sorted.calls_only

        cat polisher_output.header polisher_output.sorted.calls_only > polisher_output.vcf

        bgzip polisher_output.vcf
        tabix -p vcf polisher_output.vcf.gz

        mkdir polisher_vcf_output

        cp polisher_output.vcf.gz polisher_vcf_output/polisher_output.no_filters.vcf.gz
        echo "~{customGQFilter}" "custom GQ filter"
        echo "~{useOptimalGQFilter}" "useOptimalGQFilter"
        # if GQ filter not passed in, check if use useOptimalGQFilter is true
        if [ -z "~{customGQFilter}" ]; then
            echo "customGQ filter not set"
            if [ "~{useOptimalGQFilter}" = true ]; then
                echo "filtering with optimal GQ filter"
                bcftools view -Oz -i 'FORMAT/GQ>20 && (ILEN = 1)' polisher_output.vcf.gz > polisher_output.GQ20_INS1.vcf.gz
                tabix -p vcf polisher_output.GQ20_INS1.vcf.gz

                bcftools view -Oz -i 'FORMAT/GQ>12 && (ILEN = -1)' polisher_output.vcf.gz > polisher_output.GQ12_DEL1.vcf.gz
                tabix -p vcf polisher_output.GQ12_DEL1.vcf.gz

                bcftools view -Oz -e 'FORMAT/GQ<=5 || (ILEN = 1) || (ILEN = -1)' polisher_output.vcf.gz > polisher_output.GQ5.notINS1orDEL1.vcf.gz
                tabix -p vcf polisher_output.GQ5.notINS1orDEL1.vcf.gz

                bcftools concat -a -Oz polisher_output.GQ20_INS1.vcf.gz polisher_output.GQ12_DEL1.vcf.gz polisher_output.GQ5.notINS1orDEL1.vcf.gz > ./polisher_vcf_output/polisher_output.vcf.gz
                tabix -p vcf ./polisher_vcf_output/polisher_output.vcf.gz

                echo `zcat polisher_output.GQ20_INS1.vcf.gz | grep -v "^#" | wc -l`

            else # don't use GQ filter
                echo "customGQ filter not set, dont use a GQ filter"
                cp polisher_output.vcf.gz polisher_vcf_output/
                cp polisher_output.vcf.gz.tbi polisher_vcf_output/
            fi
        # use single GQ filter passed in
        else
            echo "use custom GQ filter passed in"
            bcftools view -Oz ~{customGQFilter} polisher_output.vcf.gz > polisher_vcf_output/polisher_output.vcf.gz
            tabix -p vcf polisher_vcf_output/polisher_output.vcf.gz
        fi

        # vcfFile=polisher_vcf_output/polisher_output.vcf.gz contains whatever filtered were passed in
        # noFiltersPolisherVcf is the original polishing output without filters

        echo `zcat polisher_vcf_output/polisher_output.vcf.gz | grep -v "^#" | wc -l`
        echo `zcat polisher_vcf_output/polisher_output.no_filters.vcf.gz | grep -v "^#" | wc -l`
        >>>

        output {
            File vcfFile="polisher_vcf_output/polisher_output.vcf.gz"
            File vcfFileTbi="polisher_vcf_output/polisher_output.vcf.gz.tbi"
            File noFiltersPolisherVcf="polisher_vcf_output/polisher_output.no_filters.vcf.gz"
        }
        runtime {
            memory: memSizeGB + " GB"
            cpu: threadCount
            disks: "local-disk " + diskSizeGB + " SSD"
            docker: dockerImage
            preemptible: 1
        }
}
