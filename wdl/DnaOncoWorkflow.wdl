#================================================================================
#
# FILE: DnaOncoWorkflow.wdl
#
# DESCRIPTION: NGeneBio DNA ONCOaccuPanel master workflow
# Note:
# REQUIREMENTS:
# BUGS:
# NOTES:
# AUTHOR: Wonseok Yoo wonseok.yoo@ngenebio.com
# COMPANY: ngenebio
# VERSION: 0.1
# CREATED: 2022.04.24
# REVISION:
# RUN COMMAND: /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/cromwell-50.jar run /opt/ngenebio/pipeline/wdl/DnaOncoWorkflow.wdl -i /opt/ngenebio/test/testOnco.json
# VALIDATION : /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/womtool-50.jar validate /opt/ngenebio/pipeline/wdl/DnaOncoWorkflow.wdl -i /opt/ngenebio/test/testOnco.json
#================================================================================

import '/opt/ngenebio/pipeline/wdl/Setting.wdl' as Setting
import '/opt/ngenebio/pipeline/wdl/DnaFastqValidation.wdl' as FastqValidation
import '/opt/ngenebio/pipeline/wdl/DnaOncoAlignment.wdl' as DnaOncoAlignment
import '/opt/ngenebio/pipeline/wdl/DnaOncoSnv.wdl' as DnaOncoSnv
import '/opt/ngenebio/pipeline/wdl/DnaOncoCnv.wdl' as DnaOncoCnv
import '/opt/ngenebio/pipeline/wdl/DnaOncoSv.wdl' as DnaOncoSv

workflow DnaOncoWorkflow {

    Array[String] sampleNames
    Array[String] sampleFastqR1
    Array[String] sampleFastqR2
    Array[Int] sampleIds
    Array[String] sampleDisease
    Array[Float] samplePurity
    
    String inputDir
    Int runId
    String runName
    String pipelineCode
    String pipelineName
    String platform
    
    # Set Parameter
    String inputFasta = '/opt/ngenebio/dependencies/hg19/human_g1k_v37.fasta'
    String runName_ = sub(runName, ' ', '_')
    String pipelineName_ = sub(pipelineName, ' ', '_')
    String platform_ = sub(platform, ' ', '_')
    String assayDir = '/opt/ngenebio/pipeline/assay_reference/ONCO/'
    
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
   String pipelineVersion = pipelineInfo[2]

    call Setting.setOncoTarget {
        input:
            outputDir = setDir.assay,
            assayDir = assayDir,
            pipelineVersion = pipelineVersion
    }

    scatter (i in range(length(sampleIds))) {

        call DnaOncoAlignment.DnaOncoAlignment { 
            input: 
                inputFasta = inputFasta,
                inputFastqR1 = DnaFastqValidation.fastqR1[i],
                inputFastqR2 = DnaFastqValidation.fastqR2[i],
                sampleName = sampleNames[i],
                sampleId = sampleIds[i],
                runId = runId,
                logDir = setDir.logs,
                resultDir = setDir.result,
                fastqDir = setDir.fastq,
                alignDir = setDir.align,
                tmpDir = setDir.tmp,
                inputIntervalBaits = setOncoTarget.baitsFile
        }

        call DnaOncoSnv.DnaOncoSnv {
            input:
                logDir = setDir.logs,
                resultDir = setDir.result,
                variantDir = setDir.snv,
                assayDir = setDir.assay,
                alignDir = setDir.align,
                tmpDir = setDir.tmp,
                inputFinalBam = DnaOncoAlignment.finalBamFile,
                sampleName = sampleNames[i],
                sampleId = sampleIds[i],
                runId = runId,
                sampleDisease = sampleDisease[i],
                inputFasta = inputFasta,
                inputTargetAllGene = setOncoTarget.snvAllBed,
                inputTargetHrd = setOncoTarget.snvHrdBed,
                inputTargetWhitelist = setOncoTarget.snvWhiBed,
                inputTargetWhitelistTxt = setOncoTarget.snvWhiTxt,
                inputSnvFormat = setOncoTarget.snvTmpVcf
        }
    
        call DnaOncoCnv.DnaOncoCnv {
            input:
                logDir = setDir.logs,
                variantDir = setDir.cnv,
                resultDir = setDir.result,
                inputFinalBam = DnaOncoAlignment.finalBamFile,
                sampleName = sampleNames[i],
                sampleId = sampleIds[i],
                runId = runId,
                sampleDisease = sampleDisease[i],
                samplePurity = samplePurity[i],
                referenceFile = setOncoTarget.cnvReference,
                hrdGeneList = setOncoTarget.hrdGeneList,
                mmrGeneList = setOncoTarget.mmrGeneList
        }
    
        call DnaOncoSv.DnaOncoSv {
            input:
                logDir = setDir.logs,
                variantDir = setDir.sv,
                resultDir = setDir.result,
                assayDir = setDir.assay,
                inputFinalBam = DnaOncoAlignment.finalBamFile,
                referenceDir = setOncoTarget.svReference,
                sampleName = sampleNames[i],
                sampleId = sampleIds[i],
                runId = runId,
                sampleDisease = sampleDisease[i],
                targetBed = setOncoTarget.svTagetBed,
                targetTxt = setOncoTarget.svTagetTxt,
                hrdGeneList = setOncoTarget.hrdGeneList,
                mmrGeneList = setOncoTarget.mmrGeneList
        }
    }

    output {
        Array[String] finalbamFile = DnaOncoAlignment.finalBamFile
        Array[String] finalSnvFile = DnaOncoSnv.finalVcfFile
        Array[String] finalCnvFile = DnaOncoCnv.finalTsvFile
        Array[String] finalSvFile = DnaOncoSv.finalTsvFile
    }
}