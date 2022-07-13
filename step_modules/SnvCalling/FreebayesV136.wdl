#================================================================================
#
# FILE: FreebayesV136.wdl
#
# DESCRIPTION: Ditect SNV/INDEL using FreeBayes (v1.3.6)
# Note:
# REQUIREMENTS:
# BUGS:
# NOTES:
# AUTHOR: Wonseok Yoo wonseok.yoo@ngenebio.com
# COMPANY: ngenebio
# VERSION: 0.1
# CREATED: 2022.06.13
# REVISION:
# RUN COMMAND: /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/cromwell-50.jar run /opt/ngenebio/pipeline/wdl/DnaBrcaSnv.wdl -i /opt/ngenebio/test/testBrcaSnv.json
# VALIDATION : /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/womtool-50.jar validate /opt/ngenebio/pipeline/wdl/DnaBrcaSnv.wdl -i /opt/ngenebio/test/testBrcaSnv.json
#================================================================================

task runFreebayes {

    # Input Parameter
    String stepName
    String logFile
    String inputFile
    String outputDir
    String outputName
    String referenceFasta
    String roiBed
    Boolean gvcf # true and false
    
    Int? minCoverage_override # penalty for a mismatch
    Int minCoverage = select_first([minCoverage_override, 20])
    Float? minAltFrac_override # gap open penalties for deletions and insertions
    Float minAltFrac = select_first([minAltFrac_override, 0.2])

    String outputFile = outputDir + outputName + '.raw.vcf'

    # Used tools
    String freebayes = '/opt/ngenebio/app/freebayes/freebayes-1.3.6'
    String gvcfOption = if gvcf then '--gvcf --gvcf-chunk -1' else ''

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
            log_progress 'Already finished to detect SNV/INDEL using FreeBayes v1.3.6'
        else
            log_progress 'Start to detect SNV/INDEL using FreeBayes v1.3.6'

            wrap ${freebayes} ${gvcfOption} \
                              --targets ${roiBed} \
                              --fasta-reference ${referenceFasta} \
                              --min-coverage ${minCoverage} \
                              --min-alternate-fraction ${minAltFrac} \
                              ${inputFile} \
                              > ${outputFile}

            wrap md5sum ${outputFile} > '${outputFile}.md5'
            log_progress 'Finish to detect SNV/INDEL using FreeBayes v1.3.6'
        fi
    >>>

    output {
        String vcfFile = outputFile
    }
}
