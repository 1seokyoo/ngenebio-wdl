#================================================================================
#
# FILE: BreaKmer.wdl
#
# DESCRIPTION: SV calling using BreaKmer
# Note:
# REQUIREMENTS:
# BUGS:
# NOTES:
# AUTHOR: Wonseok Yoo wonseok.yoo@ngenebio.com
# COMPANY: ngenebio
# VERSION: 0.1
# CREATED: 2022.06.22
# REVISION:
# RUN COMMAND: 
# VALIDATION : 
#================================================================================

task setConfigure {

    # Input Parameter
    String stepName
    String logFile
    String sampleName
    String assayDir
    String outputDir
    String targetFile
    String bamFile
    String referenceDir
    
    String sampleDir = outputDir + sampleName + '/'
    String outputFile = assayDir + sampleName + '.breakmer.conf'

    # Used tools
    String rawFile = '/opt/ngenebio/pipeline/assay_reference/ONCO/SV/breakmer.conf'

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
            log_progress 'Already finished to set Breakmer configure file'
        else
            log_progress 'Start to set Breakmer configure file'

            wrap cp ${rawFile} ${outputFile}
            sed -i "s%<analysis_name>%${sampleName}%g" ${outputFile}
            sed -i "s%<targets_bed_file>%${targetFile}%g" ${outputFile}
            sed -i "s%<sample_bam_file>%${bamFile}%g" ${outputFile}
            sed -i "s%<analysis_dir>%${sampleDir}%g" ${outputFile}
            sed -i "s%<log_file>%${logFile}%g" ${outputFile}
            sed -i "s%<reference_data_dir>%${referenceDir}%g" ${outputFile}

            wrap md5sum ${outputFile} > '${outputFile}.md5'
            log_progress 'Finish to set Breakmer configure file'
        fi
    >>>

    output {
        String confFile = outputFile
    }
}

task runBreakmer {

    # Input Parameter
    String stepName
    String logFile
    String sampleName
    String inputFile
    String outputDir
    
    String sampleDir = outputDir + sampleName + '/'
    String outputFile = sampleDir + 'output/' + sampleName + '_summary.out'
    
    # Used tools
    String script = '/opt/ngenebio/app/BreaKmer/breakmer.py'

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
            log_progress 'Already finished to call structural variation using BreaKmer'
        else
            log_progress 'Start to call structural variation using BreaKmer'
            
            wrap source /sv-pipeline/bin/activate

            wrap python ${script} ${inputFile}

            wrap deactivate

            wrap md5sum ${outputFile} > '${outputFile}.md5'
            log_progress 'Finish to call structural variation using BreaKmer'
        fi

    >>>

    output {
        String txtFile = outputFile
    }
}

task svFilter {

    # Input Parameter
    String stepName
    String logFile
    String sampleName
    String inputFile
    String outputDir
    String hrdList
    String mmrList
    String targetFile
    String panelOfNormal
    String disease

    String sampleDir = outputDir + sampleName + '/'
    String outputFile = outputDir + sampleName + '.sv.tsv'

    # Used tools
    String script = '/opt/ngenebio/pipeline/step_modules/SvCalling/onco_svFilter.py'
    
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
            log_progress 'Already finished to filter structural variation'
        else
            log_progress 'Start to filter structural variation'
            wrap python ${script} --sample_name ${sampleName} \
                                  --input_dir ${sampleDir} \
                                  --output ${outputFile} \
                                  --hrd_file ${hrdList} \
                                  --mmr_file ${mmrList} \
                                  --target_file ${targetFile} \
                                  --pon ${panelOfNormal} \
                                  --disease ${disease}

            wrap md5sum ${outputFile} > '${outputFile}.md5'
            log_progress 'Finish to filter structural variation'
        fi
    >>>

    output {
        String tsvFile = outputFile
    }
}
