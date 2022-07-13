#================================================================================
#
# FILE: FastqQc.wdl
#
# DESCRIPTION: NGeneBio fastq QC task
# Note:
# REQUIREMENTS:
# BUGS:
# NOTES:
# AUTHOR: Wonseok Yoo wonseok.yoo@ngenebio.com
# COMPANY: ngenebio
# VERSION: 1.0
# CREATED: 2022.04.24
# REVISION:
# RUN_COMMAND: /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/cromwell-50.jar run /opt/ngenebio/pipeline/step_modules/FastqValidation/FastqQc.wdl -i /opt/ngenebio/pipeline/step_modules/FastqValidation/testFastqQc.json
# VALIDATION: /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/womtool-50.jar validate /opt/ngenebio/pipeline/step_modules/FastqValidation/FastqQc.wdl
#================================================================================

workflow FastqQcRun {
    
    String readStrand 
    String logFile 
    String inputFastq 
    String inputName 
    String outputDir 
    Int cpu 

    call FastqQcTask {
        input:
            readStrand = readStrand,
            logFile = logFile,
            inputFastq = inputFastq,
            inputName = inputName,
            outputDir = outputDir,
            cpu = cpu
    }
    
}

task FastqQcTask {

    # Input Parameter
    String readStrand
    String stepName = 'FastqQC_' + readStrand
    String logFile

    String inputFastq
    String inputName
    String outputDir    
    Int cpu

    # Set Parameter
    String fastqcTool = '/opt/ngenebio/app/FastQC-0.11.3/fastqc'
    String outputFastqc = outputDir + sub(basename(inputName), '\\.fastq.gz$', '_fastqc.zip')
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

        ret=$(md5sum_check ${outputFastqc})
        if [[ $ret = 'True' ]]; then
            log_progress 'Already finished to qc fastq file (${readStrand})'
        else
            log_progress 'Start to qc fastq file (${readStrand})'

            wrap ${fastqcTool} -t ${cpu} -o ${outputDir} --extract -f fastq ${inputFastq}
            wrap md5sum ${outputFastqc} > ${outputFastqc}.md5

            log_progress 'Finished to qc fastq file (${readStrand})'
        fi
    >>>

    output {
        String fastqcFile = outputFastqc
    }
}
