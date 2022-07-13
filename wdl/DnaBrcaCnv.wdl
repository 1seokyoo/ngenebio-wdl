#================================================================================
#
# FILE: DnaBrcaCnv.wdl
#
# DESCRIPTION: NGeneBio DNA BRCAaccuTest DnaBrcaCnv workflow
# Note:
# REQUIREMENTS:
# BUGS:
# NOTES:
# AUTHOR: Wonseok Yoo wonseok.yoo@ngenebio.com
# COMPANY: ngenebio
# VERSION: 0.1
# CREATED: 2022.06.13
# REVISION:
# RUN COMMAND: /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/cromwell-50.jar run /opt/ngenebio/pipeline/wdl/DnaBrcaCnv.wdl -i /opt/ngenebio/test/testBrcaCnv.json
# VALIDATION : /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/womtool-50.jar validate /opt/ngenebio/pipeline/wdl/DnaBrcaCnv.wdl -i /opt/ngenebio/test/testBrcaCnv.json
#================================================================================

import '/opt/ngenebio/pipeline/step_modules/CnvCalling/BrcaCnv.wdl' as BrcaCnv
import '/opt/ngenebio/pipeline/wdl/Setting.wdl' as Setting

workflow DnaBrcaCnv {

    # Input Parameter
    String stepName = 'CnvCalling'
    Array[String] sampleNames
    Array[Int] sampleIds
    Array[String] sampleBams

    Int runId
    String assayDir
    String resultDir
    String variantDir
    String logDir

    String cnvReferenceDir
    
    String ampliconBed
    String primerBed
    
    String logFile = logDir + runId + '.log.cnv.stderr.txt'

    scatter (i in range(length(sampleIds))) {
        String sampleLogFile1 = logDir + sampleIds[i] + '.' + runId + '.log.cnv.stderr.txt'

        call BrcaCnv.calculateDepth {
            input:
                stepName = stepName + '-calculateDepth',
                logFile = sampleLogFile1,
                inputFile = sampleBams[i],
                outputDir = variantDir,
                outputName = sampleNames[i],
                primerBed = primerBed
        }
    }

    Array[String] depthFiles = calculateDepth.depthFile

    call BrcaCnv.setBrcaExon {
        input:
            stepName = stepName + '-setBrcaExon',
            logFile = logFile,
            logDir = logDir,
            runId = runId,
            sampleIds = sampleIds,
            ampliconBed = ampliconBed,
            outputDir = variantDir,
            assayDir = assayDir
    }

    call BrcaCnv.calculateCopyNumber {
        input:
            stepName = stepName + '-calculateCopyNumber',
            logFile = logFile,
            logDir = logDir,
            runId = runId,
            sampleIds = sampleIds,
            inputFiles = depthFiles,
            outputDir = variantDir,
            ampliconBed = ampliconBed,
            referenceDir = cnvReferenceDir
    }

    scatter (i in range(length(sampleIds))) {
        String sampleLogFile2 = logDir + sampleIds[i] + '.' + runId + '.log.cnv.stderr.txt'
        String copyNumberFile = variantDir + sampleNames[i] + '.copynumber.txt'
        String finalOutput = resultDir + sampleNames[i] + '.cnv.tsv'

        call BrcaCnv.callCnv {
            input:
                stepName = stepName + '-callCnv',
                logFile = sampleLogFile2,
                inputFile = copyNumberFile,
                tmpFile = calculateCopyNumber.txtFile,
                domainFile = setBrcaExon.txtFile,
                outputDir = variantDir,
                outputName = sampleNames[i],
                copyLossCutoff = 0.2,
                copyGainCutoff = 0.25,
                exonThrshold = 0.6
        }

        call Setting.copyFinalFile {
            input :
                stepName = stepName + '-copyFinal',
                logFile = logFile,
                inputFile = callCnv.txtFile,
                outputName = finalOutput
        }
    }

    # Outputs that will be retained when execution is complete
    output {
        Array[String] finalFile = copyFinalFile.outputFile
    }
}