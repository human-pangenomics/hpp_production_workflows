version 1.0

import "../tasks/parliament.wdl" as runParliament
import "../tasks/sniffles.wdl" as runSniffles
import "../tasks/filterSV.wdl" as runFilterSV
import "../tasks/iris.wdl" as runIris
import "../tasks/jasmine.wdl" as runJasmine

workflow sv_assembly{
    
    input{
        File ParlInputBam
        File ParlRefGenome
        File ParlIndexBam
        File ParlIndexGenome
        String? ParlPrefix
        Boolean? ParlFilterShortContigs
        Boolean? RunParl = false
        String? ParlOtherArgs
        String? ParlDockerImage
        File HifiInputBam
        File OntInputBam
        File GenomeIn
        File ReadsIn
        String HifiSnifflesOutputName
        String? SnifflesDockerImage
        String OntSnifflesOutputName
        String HifiFilterOutputName
        String OntFilterOutputName
        String? FilterDockerImage
        String HiFiIrisOutputName
        String IrisOut
        String? IrisDockerImage
        String JasmineOutputName
        String? JasmineDockerImage
        Int? maxDist
        Float? minSeqID
        Int? specReads
    }

    # Run PARLIAMENT on Illumina data if user has chosen to

    if (RunParl == true){
        call runParliament.Parliament as Parl{
            input:
                inputBam = ParlInputBam,
                refGenome = ParlRefGenome,
                indexBam = ParlIndexBam,
                indexGenome = ParlIndexGenome,
                prefix = ParlPrefix,
                filterShortContigs = ParlFilterShortContigs,
                otherArgs = ParlOtherArgs,
                dockerImage = ParlDockerImage
        }
    }

    # Run SNIFFLES on HIFI data

    call runSniffles.Sniffles as HiFiSniffles{
        input:
            inputBam = HifiInputBam,
            outputName = HifiSnifflesOutputName,
            dockerImage = SnifflesDockerImage
    }
    
    # Run SNIFFLES on ONT data

    call runSniffles.Sniffles as OntSniffles{
        input:
            inputBam = OntInputBam,
            outputName = OntSnifflesOutputName,
            dockerImage = SnifflesDockerImage
    }

    # Run filter.py on HiFi Sniffles Output

    call runFilterSV.Filter as HiFiFilter{
        input:
            inputVcf = HiFiSniffles.outputFile,
            outputName = HifiFilterOutputName,
            dockerImage = FilterDockerImage
    }

    # Run filter.py on Ont Sniffles Output

    call runFilterSV.Filter as OntFilter{
        input:
            inputVcf = OntSniffles.outputFile,
            outputName = OntFilterOutputName,
            dockerImage = FilterDockerImage
    }

    call runIris.Iris as Iris{
        input:
            genomeIn = GenomeIn,
            readsIn = ReadsIn,
            vcfIn = HiFiFilter.outputFile,
            vcfOut = HiFiIrisOutputName,
            IrisOut = IrisOut,
            dockerImage = IrisDockerImage
    }

    call runJasmine.Jasmine as Jasmine{
        input:
            InputVCFs = select_all([Iris.outputFile, OntFilter.outputFile, Parl.ParliamentVCF]),
            SV_like_errors = JasmineOutputName,
            dockerImage = JasmineDockerImage,
            maxDist = maxDist,
            minSeqID = minSeqID,
            specReads = specReads
    }

    output{
        File SnifflesHiFiOutput = HiFiSniffles.outputFile
        File SnifflesOntOutput = OntSniffles.outputFile
        File FilterHiFiOutput = HiFiFilter.outputFile
        File FilterOntOutput = OntFilter.outputFile
        File IrisHiFiOutput = Iris.outputFile
        File SV_filelist = Jasmine.SV_filelist
        File SV_like_errors = Jasmine.outputFile
    }
}
