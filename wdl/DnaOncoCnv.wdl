#================================================================================
#
# FILE: DnaOncoCnv.wdl
#
# DESCRIPTION: NGeneBio DNA ONCOaccuPanel DnaOncoCnv workflow
# Note:
# REQUIREMENTS:
# BUGS:
# NOTES:
# AUTHOR: Juyeon hong hong.juyeon@ngenebio.com
# COMPANY: ngenebio
# VERSION: 0.1
# CREATED: 2022.05.18
# REVISION:
# RUN COMMAND: /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/cromwell-50.jar run /opt/ngenebio/pipeline/wdl/DnaOncoCnv.wdl -i /opt/ngenebio/test/testOncoCnv.json
# VALIDATION : /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/womtool-50.jar validate /opt/ngenebio/pipeline/wdl/DnaOncoCnv.wdl -i /opt/ngenebio/test/testOncoCnv.json
#================================================================================

import '/opt/ngenebio/pipeline/step_modules/CnvCalling/Cnvkit083.wdl' as Cnvkit
import '/opt/ngenebio/pipeline/wdl/Setting.wdl' as Setting


workflow DnaOncoCnv {

    # Input Parameter
    String stepName = 'CnvCalling'
    String logDir
    String variantDir 
    String resultDir
    String inputFinalBam 
    String sampleName
    Int sampleId
    Int runId
    String sampleDisease
    Float samplePurity
    String referenceFile
    String hrdGeneList
    String mmrGeneList
    Int cpu = 2
    Int ram = 4

    String sampleDir = variantDir + sampleName + '/'
    String logFile = logDir + sampleId + '.' + runId + '.log.cnv.stderr.txt'
    String finalCnvOutput = resultDir + sampleName + '.cnv.tsv'

    call Setting.checkMd5Sum {
        input:
            stepName = stepName + '-checkMd5Sum',
            logFile = logFile,
            inputFile = finalCnvOutput
    }

    if (checkMd5Sum.workMd5Sum == 'true') {
        call Setting.pseudoOutput {
            input:
                inputFile = finalCnvOutput
        }
    }

    if (checkMd5Sum.workMd5Sum != 'true') {

        call Cnvkit.cnvBatch {
            input:
                stepName = stepName + '-batch',
                sampleName = sampleName,
                inputBam = inputFinalBam,
                outputDir = sampleDir,
                logFile = logFile,
                referenceFile = referenceFile,
                cpu = cpu
        }

        call Cnvkit.cnvCall {
            input:
                stepName = stepName + '-call',
                sampleName = sampleName,
                inputCns = cnvBatch.cnsFile,
                outputDir = sampleDir,
                inputPurity = samplePurity,
                logFile = logFile
        }

        call Cnvkit.cnvFilter {
            input:
                stepName = stepName + '-filter',
                sampleName = sampleName,
                inputRawCns = cnvBatch.cnsFile,
                inputCallCns = cnvCall.cnsFile,
                outputDir = sampleDir,
                inputDisease = sampleDisease,
                hrdList = hrdGeneList,
                mmrList = mmrGeneList,
                logFile = logFile
        }

        call Setting.copyFinalFile {
            input :
                stepName = stepName + '-copyFinalFile',
                logFile = logFile,
                inputFile = cnvFilter.tsvFile,
                outputName = finalCnvOutput
        }
    }

    output {
        String finalTsvFile = select_first([copyFinalFile.outputFile, pseudoOutput.pseudoFile])
    }
}