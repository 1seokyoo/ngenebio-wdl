#================================================================================
#
# FILE: DnaBrcaSnv.wdl
#
# DESCRIPTION: NGeneBio DNA BRCAaccuTest DnaBrcaSnv workflow
# Note:
# REQUIREMENTS:
# BUGS:
# NOTES:
# AUTHOR: Wonseok Yoo wonseok.yoo@ngenebio.com
# COMPANY: ngenebio
# VERSION: 0.1
# CREATED: 2022.06.13
# REVISION:
# RUN COMMAND: 
# VALIDATION : 
#================================================================================

import '/opt/ngenebio/pipeline/step_modules/SnvCalling/FreebayesV136.wdl' as RunFreebayesV136
import '/opt/ngenebio/pipeline/step_modules/SnvCalling/VcfValidate.wdl' as VcfValidation
import '/opt/ngenebio/pipeline/step_modules/SnvCalling/VcfFilter.wdl' as VcfFilter
import '/opt/ngenebio/pipeline/wdl/Setting.wdl' as Setting

workflow DnaBrcaSnv {

    # Input Parameter
    String stepName = 'SnvCalling'
    String inputBam
    String sampleName
    Int sampleId
    Int runId
    String resultDir
    String variantDir #/opt/ngenebio/output/runId/data/snv/
    String logDir #/opt/ngenebio/output/runId/logs/
    String tmpDir
    String roiBed
    String referenceFasta
    Int cpu = 2
    Int ram = 4
    
    String logFile = logDir + sampleId + '.' + runId + '.log.snv.stderr.txt'
    String finalOutput  = resultDir + sampleName + '.snv.vcf'

    call Setting.checkMd5Sum {
        input:
            stepName = stepName + '-checkMd5Sum',
            inputFile = finalOutput,
            logFile = logFile
    }

    if (checkMd5Sum.workMd5Sum == 'true') {
        call Setting.pseudoOutput {
            input:
                inputFile = finalOutput
        }
    }

    if (checkMd5Sum.workMd5Sum != 'true') {

        call RunFreebayesV136.runFreebayes {
            input:
                stepName = stepName + '-runFreebayes',
                logFile = logFile,
                inputFile = inputBam,
                outputDir = variantDir,
                outputName = sampleName,
                referenceFasta = referenceFasta,
                roiBed = roiBed,
                gvcf = false
        }

        call VcfValidation.vcfValidate  {
            input:
                stepName = stepName + '-validation',
                inputVcf = runFreebayes.vcfFile,
                inputFasta = referenceFasta,
                sampleName = sampleName,
                outputDir = variantDir,
                mode = 'all',
                tmpDir = tmpDir,
                logFile = logFile,
                cpu = cpu,
                ram = ram
        }

        call VcfFilter.checkConsecutive {
            input:
                stepName = stepName + '-checkConsecutive',
                logFile = logFile,
                inputBamFile = inputBam,
                inputVcfFile = vcfValidate.formatVcfFile,
                outputDir = variantDir,
                outputName = sampleName,
                referenceFasta = referenceFasta
        }

        call Setting.copyFinalFile {
            input :
                stepName = stepName + '-copyFinalFile',
                logFile = logFile,
                inputFile = checkConsecutive.vcfFile,
                outputName = finalOutput
        }
    }

    # Outputs that will be retained when execution is complete
    output {
        String finalFile = select_first([copyFinalFile.outputFile, pseudoOutput.pseudoFile])
    }
}