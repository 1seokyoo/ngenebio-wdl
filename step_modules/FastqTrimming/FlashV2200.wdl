#================================================================================
#
# FILE: FlashV2200.wdl
#
# DESCRIPTION: Using Flash to merge R1 & R2 fastq
# Note:
# REQUIREMENTS:
# BUGS:
# NOTES:
# AUTHOR: Wonseok Yoo wonseok.yoo@ngenebio.com
# COMPANY: ngenebio
# VERSION: 0.1
# CREATED: 2022.05.24
# REVISION:
# RUN COMMAND: 
# VALIDATION : 
#================================================================================

task runFlash {

    # Input parameter
    String logFile
    String stepName
    String outputDir
    String sampleName
    String inputFastqR1
    String inputFastqR2
    Int trimQuality
    Int trimReadLength

    String outputFastq = outputDir + sampleName + '.extendedFrags.fastq.gz'

    # Used tools
    String flash = '/opt/ngenebio/app/FLASH2/flash2'

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
            log_progress 'Already finished to merge fastq'
        else
            log_progress 'Start to merge fastq'

            wrap ${flash} -M 150 -t 2 --allow-outies ${inputFastqR1} ${inputFastqR2} -o ${sampleName} -z -d ${outputDir}

            wrap md5sum ${outputFastq} > '${outputFastq}.md5'
            log_progress 'Finish to merge fastq'
        fi
    >>>

    output {
        String fastqFile = outputFastq
    }
}