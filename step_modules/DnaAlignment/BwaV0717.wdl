#================================================================================
#
# FILE: BwaV0717.wdl
#
# DESCRIPTION: Fastq alignment to use BWA v0.7.17
# Note:
# REQUIREMENTS:
# BUGS:
# NOTES:
# AUTHOR: Wonseok Yoo wonseok.yoo@ngenebio.com
# COMPANY: ngenebio
# VERSION: 0.1
# CREATED: 2022.05.29
# REVISION:
# RUN COMMAND:
# VALIDATION :
#================================================================================

task bwaMem {

    # Input Parameter
    String stepName
    String logFile
    String inputFile
    String outputDir
    String outputName
    String tmpDir
    String referenceFasta
    Int cpu

    String outputFile = outputDir + outputName + '.raw.sam'

    # BWA Mem Option
    Int? mismatch_override # penalty for a mismatch
    Int mismatch = select_first([mismatch_override, 4])
    Int? gap_override # gap open penalties for deletions and insertions
    Int gap = select_first([gap_override, 6])

    # Used tools
    String bwa = '/opt/ngenebio/app/BWA/bwa-0.7.17/bwa'

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

        ret=$(md5sum_check ${outputFile})
        if [[ $ret = 'True' ]]; then
            log_progress 'Already finished to align using BWA v0.7.17'
        else
            log_progress 'Start to align using BWA v0.7.17'

            wrap ${bwa} mem -M -t ${cpu} \
                            -B ${mismatch} -O ${gap} \
                            -o ${outputFile} \
                            ${referenceFasta} \
                            ${inputFile}

            wrap md5sum ${outputFile} > '${outputFile}.md5'
            log_progress 'Finish to align using BWA v0.7.17'
        fi
    >>>

    output {
        String samFile = outputFile
    }
}