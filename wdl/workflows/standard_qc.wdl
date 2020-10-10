version 1.0

import "../tasks/dipcall.wdl" as dipcall_t
import "../tasks/merqury.wdl" as merqury_t
import "../tasks/meryl.wdl" as meryl_t
import "../tasks/asmgene.wdl" as asmgene_t
import "../tasks/quast.wdl" as quast_t

workflow standardQualityControl {

    input {
        Array[File] sampleReadsILM
        Array[File] maternalReadsILM
        Array[File] paternalReadsILM
        File paternalAssembly
        File maternalAssembly
        Boolean isMaleSample
        File referenceFasta
        File geneAnnotationFile
    }

    ### Dipcall ###
    call dipcall_t.dipcall as dipcall {
        input:
            assemblyFastaPat=paternalAssembly,
            assemblyFastaMat=maternalAssembly,
            referenceFasta=referenceFasta,
            isMaleSample=isMaleSample
    }

    ### Minimap2 Gene Stats ###
    call asmgene_t.asmgene as asmgenePaternal {
        input:
            assemblyFasta=paternalAssembly,
            genesFasta=geneAnnotationFile,
            referenceFasta=referenceFasta
    }
    call asmgene_t.asmgene as asmgeneMaternal {
        input:
            assemblyFasta=maternalAssembly,
            genesFasta=geneAnnotationFile,
            referenceFasta=referenceFasta
    }

    ### Quast ###
    call quast_t.quast as quastPaternal {
        input:
            assemblyFasta=paternalAssembly,
            extraArguments="--large --est-ref-size 3100000000 --no-icarus"
    }
    call quast_t.quast as quastMaternal {
        input:
            assemblyFasta=maternalAssembly,
            extraArguments="--large --est-ref-size 3100000000 --no-icarus"
    }

    ### Merqury ###
    call meryl_t.runMeryl as meryl {
        input:
            sampleReadsILM=sampleReadsILM,
            maternalReadsILM=maternalReadsILM,
            paternalReadsILM=paternalReadsILM,
            referenceFasta=referenceFasta
    }
    call merqury_t.merqury as merqury {
        input:
            assemblyFasta=paternalAssembly,
            altHapFasta=maternalAssembly,
            kmerTarball=meryl.sampleMerylDB,
            matKmerTarball=meryl.maternalHapmer,
            patKmerTarball=meryl.paternalHapmer
    }

	output {
        File dipcallVCF = dipcall.outputVCF
        File dipcallBED = dipcall.outputBED
        File dipcallFullOutput = dipcall.outputTarball
        File paternalGeneStats = asmgenePaternal.geneStats
        File maternalGeneStats = asmgeneMaternal.geneStats
        File paternalQuastResults = quastPaternal.outputTarball
        File maternalQuastResults = quastMaternal.outputTarball
        File merylHapmerImages = meryl.hapmerImages
        File merquryQV = merqury.QV
        File merquryResults = merqury.outputTarball
	}
}


