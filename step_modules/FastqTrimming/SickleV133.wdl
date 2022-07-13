#================================================================================
#
# FILE: SickleV133.wdl
#
# DESCRIPTION: Using Siokle to trim fastq
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

task runSickle {

    # Input parameter
    String logFile
    String stepName
    String outputDir
    String sampleName
    String inputFastqR1
    String inputFastqR2
    Int trimQuality
    Int trimReadLength

    # Output parameter
    String outputFastqR1 = outputDir + sampleName + '_trimmed_R1.fastq.gz'
    String outputFastqR2 = outputDir + sampleName + '_trimmed_R2.fastq.gz'
    String outputSingleFastq = outputDir + sampleName + '_trimmed_single.fastq.gz'

    # Used tools
    String sickle = '/opt/ngenebio/app/sickle-master/sickle'

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

        ret1=$(md5sum_check ${outputFastqR1})
        ret2=$(md5sum_check ${outputFastqR2})
        if [[ $ret1 = 'True' && $ret2 = 'True' ]]; then
            log_progress 'Already finished to trim fastq using sickle'
        else
            log_progress 'Start to trim fastq using sickle'

            wrap ${sickle} pe -g -f ${inputFastqR1} -r ${inputFastqR2} \
                              -l ${trimReadLength} -q ${trimQuality} -t sanger \
                              -o ${outputFastqR1} -p ${outputFastqR2} -s ${outputSingleFastq}

            wrap md5sum ${outputFastqR1} > '${outputFastqR1}.md5'
            wrap md5sum ${outputFastqR2} > '${outputFastqR2}.md5'
            log_progress 'Finish to trim fastq using sickle'
        fi
    >>>

    output {
        String fastqR1File = outputFastqR1
        String fastqR2File = outputFastqR2
    }
}