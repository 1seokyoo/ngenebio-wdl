#================================================================================
#
# FILE: BRCA_FastqPrimer_remover.wdl
#
# DESCRIPTION: NGeneBio DNA BRCAaccuPanel master workflow
# Note:
# REQUIREMENTS:
# BUGS:
# NOTES:
# AUTHOR: Sungyoon Kim sungyoon.kim@ngenebio.com
# COMPANY: ngenebio
# VERSION: 0.1
# CREATED: 2022.05.04
# REVISION:
# RUN COMMAND: /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/cromwell-50.jar run /opt/ngenebio/pipeline/wdl/DnaBrcaPreprocessing.wdl -i /opt/ngenebio/test/testBrcaPreprocessing.json
# VALIDATION : /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/womtool-50.jar validate /opt/ngenebio/pipeline/wdl/DnaBrcaPreprocessing.wdl -i /opt/ngenebio/test/testBrcaPreprocessing.json
#================================================================================

import '/opt/ngenebio/pipeline/step_modules/FastqTrimming/BrcaPreprocessing.wdl' as BrcaPreprocessing
import '/opt/ngenebio/pipeline/step_modules/FastqTrimming/SickleV133.wdl' as SickleV133
import '/opt/ngenebio/pipeline/step_modules/FastqTrimming/FlashV2200.wdl' as FlashV2200
import '/opt/ngenebio/pipeline/wdl/Setting.wdl' as Setting

workflow DnaBrcaPreprocessing {

    # Input Parameter
    String stepName = 'Preprocessing'
    String inputFastqR1 ## copied FastqR1 (output of DnaFastqValidation)
    String inputFastqR2 ## copied FastqR2 (output of DnaFastqValidation)
    String sampleName
    String sampleNameAll = sampleName + '_all'
    Int sampleId
    Int runId
    String resultDir
    String fastqDir #/opt/ngenebio/output/runId/fastq/
    String logDir #/opt/ngenebio/output/runId/logs/
    String assayDir #/opt/ngenebio/output/runId/assay_reference/
    String primerBed

    String logFile = logDir + sampleId + '.' + runId + '.log.preprocess.stderr.txt'
    String finalOutput = resultDir + sampleName + '.extendedFrags.fastq.gz'
    String finalOutputAll = resultDir + sampleNameAll + '.extendedFrags.fastq.gz'

    call BrcaPreprocessing.normalizeReadDepth {
        input:
            stepName = stepName + '-normalizeReadDepth',
            logFile = logFile,
            outputDir = fastqDir,
            sampleName = sampleName,
            inputFastqR1 = inputFastqR1,
            inputFastqR2 = inputFastqR2,
            cycle = 151,
            expectedDepth = 400000,
            targetLength = 20881
    }

    # Removed Primer Fastq Preprocessing
    call Setting.checkMd5Sum as removePrimerCheckMd5Sum {
        input:
            stepName = stepName + '-md5Check',
            logFile = logFile,
            inputFile = finalOutput,
    }

    if (removePrimerCheckMd5Sum.workMd5Sum== 'true') {
        call Setting.pseudoOutput as removePrimerPseudoOutput {
            input:
                inputFile = finalOutput
        }
    }

    if (removePrimerCheckMd5Sum.workMd5Sum!= 'true') {

        call BrcaPreprocessing.removePrimer {
            input:
                stepName = stepName + '-removePrimer',
                logFile = logFile,
                outputDir = fastqDir,
                assayDir = assayDir,
                sampleName = sampleName,
                inputFastqR1 = normalizeReadDepth.fastqR1File,
                inputFastqR2 = normalizeReadDepth.fastqR2File,
                primerBed = primerBed,
        }

        call SickleV133.runSickle as removePrimerRunSickle {
            input:
                stepName = stepName + '-sickle',
                logFile = logFile,
                outputDir = fastqDir,
                sampleName = sampleName,
                inputFastqR1 = removePrimer.fastqR1File,
                inputFastqR2 = removePrimer.fastqR2File,
                trimQuality = 30,
                trimReadLength = 100
        }

        call FlashV2200.runFlash as removePrimerRunFlash {
            input:
                stepName = stepName + '-flash',
                logFile = logFile,
                outputDir = fastqDir,
                sampleName = sampleName,
                inputFastqR1 = removePrimerRunSickle.fastqR1File,
                inputFastqR2 = removePrimerRunSickle.fastqR2File,
                trimQuality = 30,
                trimReadLength = 100
        }

        call Setting.copyFinalFile as removePrimerCopyFinalFastq {
            input :
                stepName = stepName + '-copyFinal',
                logFile = logFile,
                inputFile = removePrimerRunFlash.fastqFile,
                outputName = finalOutput
        }
    }

    # Retained Primer Fastq Preprocessing
    call Setting.checkMd5Sum as retainPrimerCheckMd5Sum{
        input:
            stepName = stepName + '-md5CheckAll',
            logFile = logFile,
            inputFile = finalOutputAll,
    }

    if (retainPrimerCheckMd5Sum.workMd5Sum== 'true') {
        call Setting.pseudoOutput as retainPrimerPseudoOutput {
            input:
                inputFile = finalOutputAll
        }
    }

    if (retainPrimerCheckMd5Sum.workMd5Sum!= 'true') {

        call SickleV133.runSickle as retainPrimerRunSickle {
            input:
                stepName = stepName + '-sickleAll',
                logFile = logFile,
                outputDir = fastqDir,
                sampleName = sampleNameAll,
                inputFastqR1 = normalizeReadDepth.fastqR1File,
                inputFastqR2 = normalizeReadDepth.fastqR2File,
                trimQuality = 30,
                trimReadLength = 100
        }

        call FlashV2200.runFlash as retainPrimerRunFlash {
            input:
                stepName = stepName + '-flashAll',
                logFile = logFile,
                outputDir = fastqDir,
                sampleName = sampleNameAll,
                inputFastqR1 = retainPrimerRunSickle.fastqR1File,
                inputFastqR2 = retainPrimerRunSickle.fastqR2File,
                trimQuality = 30,
                trimReadLength = 100
        }

        call Setting.copyFinalFile as retainPrimerCopyFinalFastq {
            input :
                stepName = stepName + '-copyFinalAll',
                logFile = logFile,
                inputFile = retainPrimerRunFlash.fastqFile,
                outputName = finalOutputAll
        }
    }

    # Outputs that will be retained when execution is complete
    output {
        String finalRemovePrimerFastqFile = select_first([removePrimerCopyFinalFastq.outputFile, removePrimerPseudoOutput.pseudoFile])
        String finalRetainPrimerFastqFile = select_first([retainPrimerCopyFinalFastq.outputFile, retainPrimerPseudoOutput.pseudoFile])
    }
}