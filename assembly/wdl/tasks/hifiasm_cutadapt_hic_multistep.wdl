version 1.0


workflow runHiCHifiasm{
    call hicHifiasm
}

# hifiasm HiC steps
# 1st: pass HiFi / not need to pass UL and HiC (extraOptions="--bin-only")
# 2nd: pass UL and fake HiFi / not need to pass HiC (extraOptions="--bin-only")
# 3rd: pass UL, HiC and fake HiFi
task hicHifiasm {
    input{
        Array[File]? childReadsHiC1
        Array[File]? childReadsHiC2
        Array[File] childReadsHiFi
        Array[File]? childReadsUL
        Int? homCov
        String childID
	String? extraOptions
        File? inputBinFilesTarGz
        # runtime configurations
        Int memSizeGB
        Int threadCount
        Int diskSizeGB
        Int preemptible
        String dockerImage
        String zones
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        ## if bin files are given we have to extract them in the directory where hifiasm is being run
        ## this enables hifiasm to skip the time-consuming process of finding read overlaps
        if [ ! -v ~{inputBinFilesTarGz} ]; then
            tar -xzf ~{inputBinFilesTarGz} --strip-components 1
            rm -rf ~{inputBinFilesTarGz} || true
        fi

        ## run trio hifiasm https://github.com/chhylp123/hifiasm
        # If ONT ultra long reads are provided
        if [[ -n "~{sep="" childReadsUL}" ]]; then
            if [[ -n "~{sep="" childReadsHiC1}" ]]; then                 
                # Keep only necessary bin files for step 3
                mkdir kept_bin_files
                mv *.ec.bin *.ovlp.reverse.bin *.ovlp.source.bin *.ul.ovlp.bin kept_bin_files
                rm -rf *.bin
                mv kept_bin_files/* .
                rm -rf kept_bin_files

                # Run step 3
                hifiasm ~{extraOptions} -o ~{childID} --ul ~{sep="," childReadsUL} --hom-cov ~{homCov} -t~{threadCount} --h1 "~{sep="," childReadsHiC1}" --h2 "~{sep="," childReadsHiC2}"  ~{sep=" " childReadsHiFi}
            else
                # Keep only necessary bin files for step 2
                mkdir kept_bin_files
                mv *.ec.bin *.ovlp.reverse.bin *.ovlp.source.bin kept_bin_files
                rm -rf *.bin
                mv kept_bin_files/* .
                rm -rf kept_bin_files

                # Run step 2
                hifiasm ~{extraOptions} -o ~{childID} --ul ~{sep="," childReadsUL} --hom-cov ~{homCov} -t~{threadCount} ~{sep=" " childReadsHiFi}
            fi
        else  
            # Run step 1
            hifiasm ~{extraOptions} -o ~{childID} -t~{threadCount} ~{sep=" " childReadsHiFi}
        fi

        #Move bin and gfa files to saparate folders and compress them 
        mkdir ~{childID}.raw_unitig_gfa
        mkdir ~{childID}.hap1.contig_gfa
        mkdir ~{childID}.hap2.contig_gfa
        mkdir ~{childID}.binFiles

        # before hardlinking gfa files to the corresponding directory make sure they exist
        # first and second step of hifiasm does not output gfa files
        if [[ -n $(find . -maxdepth 1 -iname "*.hap1.p_ctg.*") ]];
        then
            ln *.r_utg.* ~{childID}.raw_unitig_gfa
            ln *.hap1.p_ctg.* ~{childID}.hap1.contig_gfa
            ln *.hap2.p_ctg.* ~{childID}.hap2.contig_gfa
        else
            # To avoid making a separete new task for steps 1 and 2 of hifiasm
            # we make empty gfa files since we cannot have optional outputs in wdl
            touch empty.hap1.p_ctg.gfa
            touch empty.hap2.p_ctg.gfa
        fi

        ln *.bin ~{childID}.binFiles
        
        
        # make archives
        tar -cf ~{childID}.raw_unitig_gfa.tar ~{childID}.raw_unitig_gfa
        tar -cf ~{childID}.hap1.contig_gfa.tar ~{childID}.hap1.contig_gfa
        tar -cf ~{childID}.hap2.contig_gfa.tar ~{childID}.hap2.contig_gfa
        tar -cf ~{childID}.binFiles.tar ~{childID}.binFiles
        
        # compress
        pigz -p~{threadCount} ~{childID}.raw_unitig_gfa.tar
        pigz -p~{threadCount} ~{childID}.hap1.contig_gfa.tar
        pigz -p~{threadCount} ~{childID}.hap2.contig_gfa.tar
        pigz -p~{threadCount} ~{childID}.binFiles.tar
    >>>

    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        preemptible : preemptible
        cpuPlatform: "Intel Cascade Lake"
        zones : zones
    }

    output {
        File outputHaplotype1Gfa = glob("*.hap1.p_ctg.gfa")[0]
        File outputHaplotype2Gfa = glob("*.hap2.p_ctg.gfa")[0]
        File outputHaplotype1ContigGfa = "~{childID}.hap1.contig_gfa.tar.gz"
        File outputHaplotype2ContigGfa = "~{childID}.hap2.contig_gfa.tar.gz"
        File outputRawUnitigGfa = "~{childID}.raw_unitig_gfa.tar.gz"
        File outputBinFiles = "~{childID}.binFiles.tar.gz"
    }
}

