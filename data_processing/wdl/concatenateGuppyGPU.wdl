version 1.0

workflow parallelGuppyGPU {
	input {
		# input must be tar files
		Array[File] fast5_tar_files
		
	}


	# scatter in case multiple tar files are given
	scatter (fast5_tar in fast5_tar_files) {
		call splitFast5s {
			input:
				file_to_split_tar = fast5_tar

		}

		# call guppyGPU on each of the smaller "proportioned" tar files
		scatter (split_fast5_tar in splitFast5s.split_fast5_tar) {
			call guppyGPU {
				input:
					fast5_tar_file = split_fast5_tar
			}
		}

		# ?????
		call concatenateFiles as bamFile {
			input:
				files = guppyGPU.pass_bam,
				file_type = "bam"
		}

		call concatenateFiles as fastqFile {
			input:
				files = guppyGPU.pass_fastq,
				file_type = "fastq"
		}

		call concatenateFiles as summaryFile {
			input:
				files = guppyGPU.summary,
				file_type = "txt"
		}

	}

	# gather??
	output {
		Array[File] bams = bamFile.concatenatedFile
		Array[File] fastqs = fastqFile.concatenatedFile
		Array[File] summaries = summaryFile.concatenatedFile
	}
	
}

task concatenateFiles {
	input {
		Array[File] files
		String file_type
		String sample_name
		String guppy_version


		String dockerImage = "tpesout/megalodon:latest"

		# runtime
		Int preempts = 3
		Int memSizeGB = 8
		Int threadCount = 3
		Int diskSizeGB = 500
	}
	

	command {
		if [[ ${file_type} == "bam" ]]
		then
			samtools merge -o "${sample_name}_${guppy_version}.${file_type}" ${sep=" " files}

		elif [[ ${file_type} == "fastq" ]]
		then
			cat ${sep=" " files} > "${sample_name}_${guppy_version}.${file_type}"
			gzip "${sample_name}_${guppy_version}.${file_type}"

		else 
			cat ${sep=" " files} > "tmp.${file_type}"
			# remove duplicate headers
			awk 'NR==1 || !/^filename/' "tmp.${file_type}" > "${sample_name}_${guppy_version}.${file_type}"
		fi
	}

	output {
		File concatenatedFile = "${sample_name}_${guppy_version}.${file_type}*"
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
		File file_to_split_tar
		Int desired_size_GB

		String dockerImage = "jiminpark/guppy-wdl:latest" 

		# runtime
		Int preempts = 3
		Int memSizeGB = 8
		Int extraDisk = 5
		Int threadCount = 2
	}

	Int file_size = ceil(size(file_to_split_tar, "GB"))
	Int diskSizeGB = 3 * file_size + extraDisk


	command <<<
		## Extract tar file to 
		mkdir tmp
		
		# place all extracted files into directory tmp
		tar xvf "~{file_to_split_tar}" --directory tmp
		rm ~{file_to_split_tar}

		mkdir input

		# filter out only fast5 files
		for FILE in `find tmp/ -name "*.fast5"`
		do
			mv $FILE input
		done


		OUTPUT_IDX=0
		OUTPUT_DIR=fast5_tar_$OUTPUT_IDX
		mkdir $OUTPUT_DIR
		for FILE in `ls input`
		do
			size=$(du -s -BG $OUTPUT_DIR | sed 's/G.*//')
			echo $size
			if (( $size > ~{desired_size_GB} ))
			then
				tar -czvf fast5_tarball_$OUTPUT_IDX.tar.gz $OUTPUT_DIR/*
				rm -r $OUTPUT_DIR
				OUTPUT_IDX=$(($OUTPUT_IDX + 1))
				OUTPUT_DIR=fast5_tar_$OUTPUT_IDX
				mkdir $OUTPUT_DIR
			fi
			echo $(du -s -BG $OUTPUT_DIR | sed 's/G.*//')
			mv input/$FILE $OUTPUT_DIR
		done
		tar -czvf fast5_tarball_$OUTPUT_IDX.tar.gz $OUTPUT_DIR/*
		rm -r $OUTPUT_DIR

	>>>

	output {
		Array[File] split_fast5_tar = glob("*tar.gz")
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
		Int threadCount = 10
		Int extraDisk = 5
		Int gpuCount = 1
		Int maxRetries = 0
		String gpuType = "nvidia-tesla-v100"
		String nvidiaDriverVersion = "418.87.00"
		String zones = "us-west1-b"
	}

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
		File pass_bam = glob("output/pass/*.bam")[0]
		File pass_fastq = glob("output/pass/*.fastq")[0]
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