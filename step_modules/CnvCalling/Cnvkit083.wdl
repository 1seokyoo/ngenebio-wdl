#================================================================================
#
# FILE: Cnvkit083.wdl
#
# DESCRIPTION: NGeneBio DNA ONCOaccuPanel CNVkit 0.8.3 workflow
# Note:
# REQUIREMENTS:
# BUGS:
# NOTES:
# AUTHOR: Juyeon hong hong.juyeon@ngenebio.com
# COMPANY: ngenebio
# VERSION: 0.1
# CREATED: 2022.05.18
# REVISION:
# RUN COMMAND: 
# VALIDATION : /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/womtool-50.jar validate /opt/ngenebio/pipeline/step_modules/CnvCalling/Cnvkit083.wdl
#================================================================================

task cnvBatch {

    # Input Parameter
    String stepName
    String sampleName
    String inputBam
    String outputDir
    String logFile
    String referenceFile
    Int cpu 

    String outputCnr = outputDir + sampleName + '.raw.cnr'
    String outputCns = outputDir + sampleName + '.raw.cns'
    String tmpCnr = outputDir + sampleName + '.final.cnr'
    String tmpCns = outputDir + sampleName + '.final.cns'
    
    # Used tools
    String script = '/opt/ngenebio/app/CNVkit/cnvkit-0.8.3/cnvkit.py'
    
    command <<<
    
        set -e
        set -o pipefail # trace ERR through pipes
        set -o errtrace  # trace ERR through 'time command' and other functions
        set -o nounset   # set -u : exit the script if you try to use an uninitialised variable
        set -o errexit   # set -e : exit the script if any statement returns a non-true return value
        ulimit -n 2048 # Open file count

        export inputStep=${stepName}
        export logFile=${logFile}
        source /opt/ngenebio/pipeline/script/validate.sh

        ret=$(md5sum_check ${outputCns})

        if [[ $ret = 'True' ]]; then
            log_progress 'Run CNVkit batch already finished!!!'

        else
            log_progress 'Run CNVkit batch start'

            wrap source /cnv-pipeline/bin/activate

            wrap python ${script} batch ${inputBam} \
                                        --output-dir ${outputDir} \
                                        --reference ${referenceFile} \
                                        --drop-low-coverage \
                                        --processes ${cpu}

            wrap mv ${tmpCnr} ${outputCnr}
            wrap mv ${tmpCns} ${outputCns}

            wrap md5sum ${outputCns} > '${outputCns}.md5'
            log_progress 'Run CNVkit batch finished'
            wrap deactivate
        fi
    >>>

    output {
        String cnrFile = outputCnr
        String cnsFile = outputCns
    }

}


task cnvCall {

    # Input Parameter
    String stepName
    String sampleName
    String inputCns
    String outputDir
    Float inputPurity
    String logFile

    String outputCns = outputDir + sampleName + '.call.cns'

    # Used tools
    String script = '/opt/ngenebio/app/CNVkit/cnvkit-0.8.3/cnvkit.py'
    
    command <<<
    
        set -e
        set -o pipefail # trace ERR through pipes
        set -o errtrace  # trace ERR through 'time command' and other functions
        set -o nounset   # set -u : exit the script if you try to use an uninitialised variable
        set -o errexit   # set -e : exit the script if any statement returns a non-true return value
        ulimit -n 2048 # Open file count

        export inputStep=${stepName}
        export logFile=${logFile}
        source /opt/ngenebio/pipeline/script/validate.sh
                
        ret=$(md5sum_check ${outputCns})
        if [[ $ret = 'True' ]]; then
            log_progress 'Run CNVkit call already finished!!!'
        else
            log_progress 'Run CNVkit call start'
            wrap source /cnv-pipeline/bin/activate
           
            wrap python ${script} call ${inputCns} \
                                       --output ${outputCns} \
                                       --male-reference \
                                       --method clonal \
                                       --purity ${inputPurity}

            wrap md5sum ${outputCns} > '${outputCns}.md5'
            log_progress 'Run CNVkit call finished'
            wrap deactivate
        fi
    >>>

    output {
        String cnsFile = outputCns
    }
}


task cnvFilter {

    String stepName
    String sampleName
    String inputRawCns
    String inputCallCns
    String outputDir
    String inputDisease
    String hrdList
    String mmrList
    String logFile

    String outputTsv = outputDir + sampleName + '.cnv.tsv'
    String inputTargetCoverage = outputDir + sampleName + '.final.targetcoverage.cnn'

    String script = '/opt/ngenebio/pipeline/step_modules/CnvCalling/onco_cnvFilter.py'

    command <<<
    
        set -e
        set -o pipefail # trace ERR through pipes
        set -o errtrace  # trace ERR through 'time command' and other functions
        set -o nounset   # set -u : exit the script if you try to use an uninitialised variable
        set -o errexit   # set -e : exit the script if any statement returns a non-true return value
        ulimit -n 2048 # Open file count

        export inputStep=${stepName}
        export logFile=${logFile}
        source /opt/ngenebio/pipeline/script/validate.sh

        ret=$(md5sum_check ${outputTsv})
        if [[ $ret = 'True' ]]; then
            log_progress 'CNV filter already finished!!!'
        else
            log_progress 'CNV filter start'

            wrap source /cnv-pipeline/bin/activate
            wrap python ${script} -i ${inputCallCns} -j ${inputRawCns} -t ${inputTargetCoverage} \
                                  -o ${outputTsv} -s ${sampleName} -d ${inputDisease} \
                                  -r ${hrdList} -m ${mmrList}

            wrap md5sum ${outputTsv} > '${outputTsv}.md5'
            log_progress 'CNV filter finished'

            deactivate
        fi

    >>>

    output {
        String tsvFile = outputTsv
    }
}