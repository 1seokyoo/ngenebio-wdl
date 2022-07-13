#================================================================================
#
# FILE: DnaBrcaWorkflow.wdl
#
# DESCRIPTION: NGeneBio DNA BRCAaccuTest master workflow
# Note:
# REQUIREMENTS:
# BUGS:
# NOTES:
# AUTHOR: Wonseok Yoo wonseok.yoo@ngenebio.com
# COMPANY: ngenebio
# VERSION: 0.1
# CREATED: 2022.05.19
# REVISION:
# RUN COMMAND: /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/cromwell-50.jar run /opt/ngenebio/pipeline/wdl/DnaBrcaWorkflow.wdl -i /opt/ngenebio/test/testBrca.json
# VALIDATION : /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/womtool-50.jar validate /opt/ngenebio/pipeline/wdl/DnaBrcaWorkflow.wdl -i /opt/ngenebio/test/testBrca.json
#================================================================================

import '/opt/ngenebio/pipeline/wdl/Setting.wdl' as Setting
import '/opt/ngenebio/pipeline/wdl/DnaFastqValidation.wdl' as FastqValidation
import '/opt/ngenebio/pipeline/wdl/DnaBrcaPreprocessing.wdl' as Preprocessing
import '/opt/ngenebio/pipeline/wdl/DnaBrcaAlignment.wdl' as Alignment
import '/opt/ngenebio/pipeline/wdl/DnaBrcaSnv.wdl' as Snv
import '/opt/ngenebio/pipeline/wdl/DnaBrcaCnv.wdl' as Cnv

workflow DnaBrcaWorkflow {

    Array[String] sampleNames
    Array[String] sampleFastqR1
    Array[String] sampleFastqR2
    Array[Int] sampleIds

    String inputDir
    Int runId
    String runName
    String pipelineCode
    String pipelineName
    String platform
    
    # Set Parameter
    String inputFasta = '/opt/ngenebio/dependencies/hg19/human_hg19.brca.fasta'
    String runName_ = sub(runName, ' ', '_')
    String pipelineName_ = sub(pipelineName, ' ', '_')
    String platform_ = sub(platform, ' ', '_')
    
    call Setting.setDir {
        input:
            runId = runId
    }

    call FastqValidation.DnaFastqValidation { 
        input: 
            inputDir = inputDir,
            runDir = setDir.run,
            fastqDir = setDir.fastq,
            fastqcDir = setDir.fastqc,
            sampleFastqR1 = sampleFastqR1,
            sampleFastqR2 = sampleFastqR2,
            sampleIds = sampleIds,
            runId = runId
    }

    call Setting.splitPipelineCode {
        input:
            pipelineCode = pipelineCode
    }
    
    Array[String] pipelineInfo = splitPipelineCode.pipelineInfo
    String panelVersion = pipelineInfo[2]
    String panelCustom = pipelineInfo[7]
    
    if (panelCustom == '0')  { String setVersionOrigin = panelVersion + '/' }
    if (panelCustom != '0')  { String setVersionCustom = panelVersion + panelCustom + '/' }
    String setVersion = select_first([setVersionCustom, setVersionOrigin])
    String assayTargetDir = '/opt/ngenebio/pipeline/assay_reference/BRCA/Target/' + setVersion
    String assayCnvDir = '/opt/ngenebio/pipeline/assay_reference/BRCA/Cnv/' + setVersion

    call Setting.setBrcaTarget {
        input:
            outputDir = setDir.assay,
            assayTargetDir = assayTargetDir
    }

    scatter (i in range(length(sampleIds))) {

        call Preprocessing.DnaBrcaPreprocessing {
            input:
                inputFastqR1 = DnaFastqValidation.fastqR1[i],
                inputFastqR2 = DnaFastqValidation.fastqR2[i],
                sampleName = sampleNames[i],
                sampleId = sampleIds[i],
                runId = runId,
                resultDir = setDir.result,
                fastqDir = setDir.fastq,
                logDir = setDir.logs,
                assayDir = setDir.assay,
                primerBed = setBrcaTarget.primerBed
        }

        call Alignment.DnaBrcaAlignment {
            input:
                inputFastq = DnaBrcaPreprocessing.finalRemovePrimerFastqFile,
                sampleName = sampleNames[i],
                sampleId = sampleIds[i],
                runId = runId,
                resultDir = setDir.result,
                alignDir  = setDir.align,
                logDir = setDir.logs,
                tmpDir = setDir.tmp,
                primerBed = setBrcaTarget.primerBed,
                primerBedpe = setBrcaTarget.primerBedpe,
                referenceFasta = inputFasta
        }

        call Snv.DnaBrcaSnv {
            input:
                inputBam = DnaBrcaAlignment.finalFile,
                sampleName = sampleNames[i],
                sampleId = sampleIds[i],
                runId = runId,
                resultDir = setDir.result,
                variantDir = setDir.snv,
                logDir = setDir.logs,
                tmpDir = setDir.tmp,
                roiBed = setBrcaTarget.roiBed,
                referenceFasta = inputFasta
        }
    }

    Array[String] cnvInputBamFiles = DnaBrcaAlignment.bqsrFile

    call Cnv.DnaBrcaCnv {
        input:
            sampleNames = sampleNames,
            sampleIds = sampleIds,
            sampleBams = cnvInputBamFiles,
            runId = runId,
            assayDir = setDir.assay,
            resultDir = setDir.result,
            variantDir = setDir.cnv,
            logDir = setDir.logs,
            cnvReferenceDir = assayCnvDir,
            ampliconBed = setBrcaTarget.ampliconBed,
            primerBed = setBrcaTarget.primerBed
    }

    output {
        Array[String] finalbamFile = DnaBrcaAlignment.finalFile
        Array[String] finalSnvFile = DnaBrcaSnv.finalFile
        Array[String] finalCnvFile = DnaBrcaCnv.finalFile
    }
}
