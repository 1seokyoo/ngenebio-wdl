#================================================================================
#
# FILE: DnaOncoSv.wdl
#
# DESCRIPTION: NGeneBio DNA ONCOaccuPanel DnaOncoSv workflow
# Note:
# REQUIREMENTS:
# BUGS:
# NOTES:
# AUTHOR: Wonseok Yoo wonseok.yoo@ngenebio.com
# COMPANY: ngenebio
# VERSION: 0.1
# CREATED: 2022.06.28
# REVISION:
# RUN COMMAND: /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/cromwell-50.jar run /opt/ngenebio/pipeline/wdl/DnaOncoSv.wdl -i /opt/ngenebio/test/testOncoSv.json
# VALIDATION : /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/womtool-50.jar validate /opt/ngenebio/pipeline/wdl/DnaOncoSv.wdl -i /opt/ngenebio/test/testOncoSv.json
#================================================================================

import '/opt/ngenebio/pipeline/step_modules/SvCalling/BreaKmer.wdl' as BreaKmer
import '/opt/ngenebio/pipeline/wdl/Setting.wdl' as Setting


workflow DnaOncoSv {

    # Input Parameter
    String stepName = 'SvCalling'
    String logDir
    String variantDir 
    String resultDir
    String assayDir
    String inputFinalBam 
    String referenceDir
    String sampleName
    Int sampleId
    Int runId
    String sampleDisease
    String targetBed
    String targetTxt
    String hrdGeneList
    String mmrGeneList

    String logFile = logDir + sampleId + '.' + runId + '.log.sv.stderr.txt'
    String finalSvOutput = resultDir + sampleName + '_sv.tsv'

    String panelOfNormal = '/opt/ngenebio/pipeline/step_modules/SvCalling/onco_pon20180201.txt'

    call Setting.checkMd5Sum {
            input:
                stepName = stepName + '-checkMd5Sum',
                logFile = logFile,
                inputFile = finalSvOutput
    }

    if (checkMd5Sum.workMd5Sum == 'true') {
        call Setting.pseudoOutput {
            input:
                inputFile = finalSvOutput
        }
    }

    if (checkMd5Sum.workMd5Sum != 'true') {

        call BreaKmer.setConfigure {
            input:
                stepName = stepName + '-setConfigure',
                logFile = logFile,
                sampleName = sampleName,
                assayDir = assayDir,
                outputDir = variantDir,
                targetFile = targetBed,
                bamFile = inputFinalBam,
                referenceDir = referenceDir
        }

        call BreaKmer.runBreakmer {
            input:
                stepName = stepName + '-runBreakmer',
                logFile = logFile,
                sampleName = sampleName,
                inputFile = setConfigure.confFile,
                outputDir = variantDir
        }

        call BreaKmer.svFilter {
            input :
                stepName = stepName + '-Filter',
                logFile = logFile,
                sampleName = sampleName,
                inputFile = runBreakmer.txtFile,
                outputDir = variantDir,
                hrdList = hrdGeneList,
                mmrList = mmrGeneList,
                targetFile = targetTxt,
                panelOfNormal = panelOfNormal,
                disease = sampleDisease
        }
        
        call Setting.copyFinalFile {
            input :
                stepName = stepName + '-copyFinalFile',
                logFile = logFile,
                inputFile = svFilter.tsvFile,
                outputName = finalSvOutput
        }
    }
            
    output {
        String finalTsvFile = select_first([copyFinalFile.outputFile, pseudoOutput.pseudoFile])
    }
        
}