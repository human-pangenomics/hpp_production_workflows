
version 1.0

workflow findHomozygousRegions {
    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "detect homozygous regions from paf alignment of two assemblies"
}
    call FindHomozygousRegions

    output{
        File bed=FindHomozygousRegions.bed
        File extendedBed=FindHomozygousRegions.extendedBed
    }
}

task FindHomozygousRegions{
    input {
        File pafFile
        String minWindowSizeBp
        String extendBp
        String outPrefix

        Int memSizeGB = 8
        Int threadCount = 2
        Int diskSizeGB = 8
        String dockerImage = "mobinasri/secphase:dev-v0.2.0"
    }
    command <<<
        set -eux -o pipefail
        set -o xtrace

        python3 /home/programs/src/find_homozygous_regions.py -p ~{pafFile} -m ~{minWindowSizeBp} -e ~{extendBp} -o ~{outPrefix}

        cut -f 4-6 ~{outPrefix}.bed | cat ~{outPrefix}.bed - | cut -f1-3 > tmp; mv tmp ~{outPrefix}.bed
        cut -f 4-6 ~{outPrefix}flanking_~{extendBp}.bed | cat ~{outPrefix}flanking_~{extendBp}.bed  - | cut -f1-3 > tmp; mv tmp ~{outPrefix}flanking_~{extendBp}.bed
  	>>>
  	output {
  		  File bed = "~{outPrefix}.bed"
  		  File extendedBed = glob("~{outPrefix}*flanking*.bed")[0]
  	}
      runtime {
          memory: memSizeGB + " GB"
          cpu: threadCount
          disks: "local-disk " + diskSizeGB + " SSD"
          docker: dockerImage
          preemptible: 1
      }
}
