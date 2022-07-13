#================================================================================
#
# FILE: BrcaCnv.wdl
#
# DESCRIPTION: BRCA CNV calling using in-house script
# Note:
# REQUIREMENTS:
# BUGS:
# NOTES:
# AUTHOR: Wonseok Yoo wonseok.yoo@ngenebio.com
# COMPANY: ngenebio
# VERSION: 0.1
# CREATED: 2022.06.14
# REVISION:
# RUN COMMAND: 
# VALIDATION : 
#================================================================================

task calculateDepth {

    # Input Parameter
    String stepName
    String logFile
    String inputFile
    String outputDir
    String outputName
    String primerBed
    
    String outputFile = outputDir + outputName + '.amplicon.depth.txt'
    String script = '/opt/ngenebio/pipeline/step_modules/CnvCalling/brca_ampliconDepth.py'

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
            log_progress 'Already finished to calculate amplicon depth using in-house script'
        else
            log_progress 'Start to calculate amplicon depth using in-house script'

            wrap python ${script} ${primerBed} ${inputFile} ${outputFile}

            wrap md5sum ${outputFile} > '${outputFile}.md5'
            log_progress 'Finish to calculate amplicon depth using in-house script'
        fi
    >>>

    output {
        String depthFile = outputFile
    }
}

task setBrcaExon {

    # Input Parameter
    String stepName
    String logFile
    String logDir
    Int runId
    Array[Int] sampleIds
    String ampliconBed
    String outputDir
    String assayDir
    
    String outputFile = assayDir + 'exon.txt'

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
            log_progress 'Already finished to parse amplicon bed making exon domain file'
        else
            log_progress 'Start to parse amplicon bed making exon domain file'

            awk '{print $6 "\t" $5 "\t" $8}' ${ampliconBed} \
             | sort -u \
             | sed -e "s/exon//g" -e "s/-$//" \
            > ${outputFile}

            wrap md5sum ${outputFile} > '${outputFile}.md5'
            log_progress 'Finish to parse amplicon bed making exon domain file'
        fi
        
        for sampleId in ${sep=' ' sampleIds}
        do
            sampleLogFile=`echo ${logDir}$sampleId'.'${runId}'.log.stderr.txt'`
            cat $logFile >> $sampleLogFile
        done
        rm $logFile
    >>>

    output {
        String txtFile = outputFile
    }
}

task calculateCopyNumber {

    # Input Parameter
    String stepName
    String logFile
    String logDir
    Int runId
    Array[Int] sampleIds
    Array[String] inputFiles
    String outputDir
    String ampliconBed
    String referenceDir

    String outputFile = outputDir + runId + '.copynumber.txt'

    # Used tools
    String script = '/opt/ngenebio/pipeline/step_modules/CnvCalling/brca_calculateCopynumber.R'

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
            log_progress 'Already finished to calculate copynumber using in-house script'
        else
            log_progress 'Start to calculate copynumber using in-house script'
            wrap Rscript ${script} ${outputDir} ${referenceDir} ${ampliconBed} ${outputFile} 

            wrap md5sum ${outputFile} > '${outputFile}.md5'
            log_progress 'Finish to calculate copynumber using in-house script'
        fi
        
        for sampleId in ${sep=' ' sampleIds}
        do
            sampleLogFile=`echo ${logDir}$sampleId'.'${runId}'.log.stderr.txt'`
            cat $logFile >> $sampleLogFile
        done
        rm $logFile
    >>>

    output {
        String txtFile = outputFile
    }
}

task callCnv {

    # Input Parameter
    String stepName
    String logFile
    String inputFile
    String tmpFile
    String domainFile
    String outputDir
    String outputName
    Float copyLossCutoff
    Float copyGainCutoff
    Float exonThrshold

    String outputFile = outputDir + outputName + '.cnv.tsv'

    # Used tools
    String script = '/opt/ngenebio/pipeline/step_modules/CnvCalling/brca_predictionCnv.py'

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
            log_progress 'Already finished to call cnv using in-house script'
        else
            log_progress 'Start to call cnv using in-house script'
            wrap python ${script} --input ${inputFile} \
                                  --output ${outputFile} \
                                  --domain ${domainFile} \
                                  --loss ${copyLossCutoff} \
                                  --gain ${copyGainCutoff} \
                                  --exon ${exonThrshold} 

            wrap md5sum ${outputFile} > '${outputFile}.md5'
            log_progress 'Finish to call cnv using in-house script'
        fi
    >>>

    output {
        String txtFile = outputFile
    }
}
