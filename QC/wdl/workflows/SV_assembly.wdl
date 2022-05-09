version 1.0

import "../tasks/parliament.wdl" as runParliament
import "../tasks/sniffles.wdl" as runSniffles
import "../tasks/filterSV.wdl" as runFilterSV
import "../tasks/iris.wdl" as runIris
import "../tasks/jasmine.wdl" as runJasmine

workflow sv_assembly{
    
    meta{
        author: "Avani Khadilkar"
        email: "akhadilk@ucsc.edu"
        description: "WDL implementation of the small variant like errors calling section of the [T2T Polishing Case Study](https://github.com/arangrhie/T2T-Polish/blob/master/doc/T2T_polishing_case_study.md)."
    }
    
    input{
        File? IlluminaInputBam
        File? IlluminaIndexBam
        File HifiInputBam
        File OntInputBam
        File RefGenome
        File IndexGenome
        String SampleName
        
        Boolean? ParlFilterShortContigs = true
        Boolean? RunParl = false
        String? ParlOtherArgs
        Int? maxDist
        Float? minSeqID
        Int? specReads
        
        String? ParlDockerImage
        String? SnifflesDockerImage
        String? FilterDockerImage
        String? IrisDockerImage
        String? JasmineDockerImage
    }

    # Run PARLIAMENT on Illumina data if user has chosen to

    if (RunParl == true){
        call runParliament.Parliament as Parl{
            input:
                inputBam = IlluminaInputBam,
                refGenome = RefGenome,
                indexBam = IlluminaIndexBam,
                indexGenome = IndexGenome,
                SampleName = SampleName,
                filterShortContigs = ParlFilterShortContigs,
                otherArgs = ParlOtherArgs,
                dockerImage = ParlDockerImage
        }
    }

    # Run SNIFFLES on HIFI data

    call runSniffles.Sniffles as HiFiSniffles{
        input:
            inputBam = HifiInputBam
    }
    
    # Run SNIFFLES on ONT data

    call runSniffles.Sniffles as OntSniffles{
        input:
            inputBam = OntInputBam
    }

    # Run filter.py on HiFi Sniffles Output

    call runFilterSV.Filter as HiFiFilter{
        input:
            inputVcf = HiFiSniffles.outputFile
    }

    # Run filter.py on Ont Sniffles Output

    call runFilterSV.Filter as OntFilter{
        input:
            inputVcf = OntSniffles.outputFile
    }

    call runIris.Iris as Iris{
        input:
            genomeIn = RefGenome,
            readsIn = HifiInputBam,
            vcfIn = HiFiFilter.outputFile
    }

    call runJasmine.Jasmine as Jasmine{
        input:
            InputVCFs = select_all([Iris.outputFile, OntFilter.outputFile, Parl.ParliamentVCF]),
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
