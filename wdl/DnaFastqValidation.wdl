
#================================================================================
#
# FILE: DnaFastqValidation.wdl
#
# DESCRIPTION: NGeneBio DNA fastq validation pipeline
# Note:
# REQUIREMENTS:
# BUGS:
# NOTES:
# AUTHOR: Wonseok Yoo wonseok.yoo@ngenebio.com
# COMPANY: ngenebio
# VERSION: 0.1
# CREATED: 2022.04.24
# REVISION:
# VALIDATION: /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/womtool-50.jar validate /opt/ngenebio/pipeline/wdl/DnaFastqValidation.wdl
#================================================================================

import '/opt/ngenebio/pipeline/step_modules/FastqValidation/FastqQc.wdl' as FastqQc
import '/opt/ngenebio/pipeline/step_modules/FastqValidation/FastqValidation.wdl' as FastqValidation
import '/opt/ngenebio/pipeline/step_modules/FastqValidation/FastqCopy.wdl' as FastqCopy
import '/opt/ngenebio/pipeline/wdl/Setting.wdl' as Setting

workflow DnaFastqValidation {

    # Input Parameter
    String stepName = 'FastqValidation'
    String inputDir
    String runDir
    String fastqDir
    String fastqcDir
    Array[String] sampleFastqR1
    Array[String] sampleFastqR2
    Array[String] sampleIds
    String runId
    Int cpu = 2

    scatter (i in range(length(sampleIds))) {

        String logFile = runDir + '/' + 'logs/' + sampleIds[i] + '.' + runId + '.log.fastqValidation.stderr.txt'

        # FastqR1
        call FastqCopy.FastqCopyTask as setFastqR1 {
            input:
                readStrand = 'R1',
                logFile = logFile,
                inputDir = inputDir,
                outputDir = fastqDir,
                inputFastq = sampleFastqR1[i]
        }

        call FastqValidation.FastqValidationTask as runFastqR1Validation {
            input:
                readStrand = 'R1',
                logFile = logFile,
                inputFastq = setFastqR1.fastqFile
        }

        call FastqQc.FastqQcTask as runFastqcR1 {
            input:
                readStrand = 'R1',
                logFile = logFile,
                inputFastq = setFastqR1.fastqFile,
                inputName = sampleFastqR1[i],
                outputDir = fastqcDir,
                cpu = cpu
        }

        # FastqR2
        call FastqCopy.FastqCopyTask as setFastqR2 {
            input:
                readStrand = 'R2',
                logFile = logFile,
                inputDir = inputDir,
                outputDir = fastqDir,
                inputFastq = sampleFastqR2[i]
        }

        call FastqValidation.FastqValidationTask as runFastqR2Validation {
            input:
                readStrand = 'R2',
                logFile = logFile,
                inputFastq = setFastqR2.fastqFile
        }

        call FastqQc.FastqQcTask as runFastqcR2 {
            input:
                readStrand = 'R2',
                logFile = logFile,
                inputFastq = setFastqR2.fastqFile,
                inputName = sampleFastqR2[i],
                outputDir = fastqcDir,
                cpu = cpu
        }

    }

    # Outputs that will be retained when execution is complete
    output {
        Array[String?] fastqR1 = setFastqR1.fastqFile
        Array[String?] fastqR2 = setFastqR2.fastqFile
        Array[String?] fastqcR1 = runFastqcR1.fastqcFile
        Array[String?] fastqcR2 = runFastqcR2.fastqcFile
    }
}