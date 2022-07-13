#================================================================================
#
# FILE: FastqCopy.wdl
#
# DESCRIPTION: NGeneBio fastq copy task
# Note:
# REQUIREMENTS:
# BUGS:
# NOTES:
# AUTHOR: Wonseok Yoo wonseok.yoo@ngenebio.com
# COMPANY: ngenebio
# VERSION: 1.0
# CREATED: 2022.04.24
# REVISION:
# RUN_COMMAND: /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/cromwell-50.jar run /opt/ngenebio/pipeline/step_modules/FastqValidation/FastqCopy.wdl -i /opt/ngenebio/pipeline/step_modules/FastqValidation/testFastqCopy.json
# VALIDATION: /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/womtool-50.jar validate /opt/ngenebio/pipeline/step_modules/FastqValidation/FastqCopy.wdl
#================================================================================
workflow FastqCopyRun {

    String readStrand
    String logFile
    String inputDir
    String outputDir
    String inputFastq

    call FastqCopyTask {
        input:
            readStrand = readStrand,
            logFile = logFile,
            inputDir = inputDir,
            outputDir = outputDir,
            inputFastq = inputFastq
    }
}


task FastqCopyTask {

    # Input Parameter
    String readStrand
    String stepName = 'FastqCopy_' + readStrand
    String logFile

    String inputDir
    String outputDir
    String inputFastq

    # Set Parameter
    String outputFastq = outputDir + inputFastq

    command <<<
    
        # Ensure script halts when non-zero exit status detected
        set -e
        set -o pipefail # trace ERR through pipes
        set -o errtrace  # trace ERR through 'time command' and other functions
        set -o nounset   # set -u : exit the script if you try to use an uninitialised variable
        set -o errexit   # set -e : exit the script if any statement returns a non-true return value
        ulimit -n 2048 # Open file count

        export inputStep=${stepName}
        export logFile=${logFile}
        source /opt/ngenebio/pipeline/script/validate.sh

        ret=$(md5sum_check ${outputFastq})
        if [[ $ret = 'True' ]]; then
            log_progress 'Already finished to copy raw file (${readStrand}) to fastq dir'
        else
            log_progress 'Start to copy raw file (${readStrand}) to fastq dir'

            wrap cp ${inputDir}/${inputFastq} ${outputFastq}
            wrap md5sum ${outputFastq} > ${outputFastq}.md5

            log_progress 'Finished to copy raw file (${readStrand}) to fastq dir'
        fi
    >>>

    output {
        String fastqFile = outputFastq
    }
}