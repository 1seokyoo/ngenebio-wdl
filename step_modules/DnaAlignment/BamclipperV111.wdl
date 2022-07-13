#================================================================================
#
# FILE: BamclipperV111.wdl
#
# DESCRIPTION: Using bamclipper to check soft-clip
# Note:
# REQUIREMENTS:
# BUGS:
# NOTES:
# AUTHOR: Wonseok Yoo wonseok.yoo@ngenebio.com
# COMPANY: ngenebio
# VERSION: 0.1
# CREATED: 2022.05.24
# REVISION:
# RUN_COMMAND: 
# VALIDATION: 
#================================================================================

task bamclipper {

    # Input Parameter
    String stepName
    String logFile
    String sampleName
    String inputFile
    String outputDir
    String primerBedpe
    Int cpu

    String outputFile = outputDir + sampleName + '.bamclip.bam'
    String bamclipOutputFile = sub(inputFile, "\\.bam", "\\.primerclipped.bam")

    # Used tools
    String bamclipper = '/opt/ngenebio/app/bamclipper/bamclipper.sh'
    String samtools = '/opt/ngenebio/app/samtools/samtools-1.3.1/samtools'
    String parallel = '/opt/ngenebio/app/parallel-20200522/src/parallel' 

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
            log_progress 'Already finished to check soft-clip using bamclipper'
        else
            log_progress 'Start to check soft-clip using bamclipper'

            wrap ${bamclipper} -b ${inputFile} -p ${primerBedpe} \
                               -s ${samtools} -g ${parallel} -n ${cpu}

            wrap mv ${bamclipOutputFile} ${outputFile}
            wrap mv ${bamclipOutputFile}.bai ${outputFile}.bai

            wrap md5sum ${outputFile} > '${outputFile}.md5'
            log_progress 'Finish to check soft-clip using bamclipper'
        fi
    >>>

    output {
        String bamFile = outputFile
    } 
}