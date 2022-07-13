#================================================================================
#
# FILE: FastqValidation.wdl
#
# DESCRIPTION: NGeneBio fastq validation task
# Note:
# REQUIREMENTS:
# BUGS:
# NOTES:
# AUTHOR: Wonseok Yoo wonseok.yoo@ngenebio.com
# COMPANY: ngenebio
# VERSION: 1.0
# CREATED: 2022.04.24
# REVISION:
# RUN_COMMAND: /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/cromwell-50.jar run /opt/ngenebio/pipeline/step_modules/FastqValidation/FastqValidation.wdl -i /opt/ngenebio/pipeline/step_modules/FastqValidation/testFastqValidation.json
# VALIDATION: /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/womtool-50.jar validate /opt/ngenebio/pipeline/step_modules/FastqValidation/FastqValidation.wdl
#================================================================================
workflow FastqValidationTaskRun {

    # Input Parameter
    String readStrand
    String logFile
    String inputFastq
    
    call FastqValidationTask {
        input:
            readStrand = readStrand,
            logFile = logFile,
            inputFastq = inputFastq
    }

}


task FastqValidationTask {

    # Input Parameter
    String readStrand
    String stepName = 'FastqValidation_' + readStrand
    String logFile

    String inputFastq

    # Set Parameter
    String fastqvalidatorTool = '/opt/ngenebio/app/fastQValidator/fastQValidator'
    String tmpFile = sub(inputFastq, '\\.fastq.gz$', '.fastq')
    String outputTxt = sub(inputFastq, '\\.fastq.gz$', '.fastq.validation.txt')

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

        ret=$(md5sum_check ${outputTxt})
        if [[ $ret = 'True' ]]; then
            log_progress 'Already finished to validate fastq file (${readStrand})'
        else
            log_progress 'Start to validate fastq file (${readStrand})'
            
            wrap gzip -dc ${inputFastq} > ${tmpFile}
            wrap ${fastqvalidatorTool} --file ${tmpFile} --minReadLen 1 --printableErrors -1 --baseSpace > ${outputTxt}
            wrap md5sum ${outputTxt} > ${outputTxt}.md5
            wrap rm -r -f ${tmpFile}

            log_progress 'Finished to validate fastq file (${readStrand})'
        fi

        if grep -F 'FASTQ_SUCCESS' ${outputTxt}; then
            echo 'TRUE' > /dev/stdout
        else
            log_error 'Failed fastq validation. You need to check the fastq file (${readStrand}).'
            echo 'FALSE' > /dev/stdout
        fi
    >>>

    output {
        # String outputTxt = read_string(stdout())
        String outputStd = stdout()
    }
}