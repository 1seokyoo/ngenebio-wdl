#================================================================================
#
# FILE: VcfFilter.wdl
#
# DESCRIPTION: NGeneBio DNA ONCOaccuPanel Vcf Filter workflow
# Note:
# REQUIREMENTS:
# BUGS:
# NOTES:
# AUTHOR: Juyeon hong hong.juyeon@ngenebio.com
# COMPANY: ngenebio
# VERSION: 0.1
# CREATED: 2022.05.10
# REVISION:
# RUN COMMAND: 
# VALIDATION : /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/womtool-50.jar validate /opt/ngenebio/pipeline/step_modules/SnvCalling/VcfFilter.wdl
#================================================================================


task vcfFilter {

    # Input Parameter
    String stepName
    String inputVcf
    String sampleName
    String databaseList  
    String vcfFormat
    String cutoff
    String clinvarFilter
    String logFile

    String outputVcf 
    String tmpVcf = sub(outputVcf , '\\.vcf$', '_tmp.vcf')
    String tmpLog = sub(outputVcf , '\\.vcf$', '_tmp.log')

    # Used tools
    String script = '/opt/ngenebio/pipeline/step_modules/SnvCalling/onco_snvFilter.py'
    String vcfAnno = '/opt/ngenebio/app/vcfanno/vcfanno'
    String vcfAnnoCallingConf = '/opt/ngenebio/pipeline/assay_reference/ONCO/COMMON/vcfanno.conf'
    String snvRecover = '/opt/ngenebio//pipeline/assay_reference/ONCO/COMMON/recoverVariant.json'

    command <<<

        # set -e
        # set -o pipefail # trace ERR through pipes
        # set -o errtrace  # trace ERR through 'time command' and other functions
        # set -o nounset   # set -u : exit the script if you try to use an uninitialised variable
        # set -o errexit   # set -e : exit the script if any statement returns a non-true return value
        # ulimit -n 2048 # Open file count

        export inputStep=${stepName}
        export logFile=${logFile}
        source /opt/ngenebio/pipeline/script/validate.sh

        ret1=$(md5sum_check ${outputVcf})
        ret2=$(vcf_check ${inputVcf})

        if [[ $ret1 = 'True' ]]; then
            log_progress 'Variant Filter already finished!!!'
        elif [[ $ret2 = 'False' ]]; then
            log_progress 'Variant Filter input vcf is empty!!!'

            wrap > ${outputVcf}
            wrap md5sum ${outputVcf} > '${outputVcf}.md5'
            log_progress 'Variant Filter finished'
        else
            log_progress 'Variant Filter start'

            # Filtering Database Annotation
            wrap source /snv-pipeline/bin/activate
            wrap ${vcfAnno} ${vcfAnnoCallingConf} ${inputVcf} > ${tmpVcf} 2> ${tmpLog}
            grep -v 'not found in INFO' ${tmpLog}

            # Variant Filter
            wrap python ${script} --input ${tmpVcf} --output ${outputVcf} --database ${databaseList} --format ${vcfFormat} --cutoff ${cutoff} --clinvar ${clinvarFilter} --recover ${snvRecover}

            wrap md5sum ${outputVcf} > '${outputVcf}.md5'
            log_progress 'Variant Filter finished'

            deactivate
        fi
    >>>

    output {
        String finalVcfFile = outputVcf
    }
}

task checkConsecutive {

    # Input Parameter
    String stepName
    String logFile
    String inputBamFile
    String inputVcfFile
    String outputDir
    String outputName
    String referenceFasta
    
    Float? minAltFrac_override # gap open penalties for deletions and insertions
    Float minAltFrac = select_first([minAltFrac_override, 0.2])

    String outputFile = outputDir + outputName + '.consecutive.vcf'

    # Used tools
    String script = '/opt/ngenebio/pipeline/step_modules/SnvCalling/brca_checkConsecutive.py'
    String consecutiveDb = '/opt/ngenebio/pipeline/assay_reference/BRCA/consecutive_variant.txt'
    String freebayes = '/opt/ngenebio/app/freebayes/freebayes-1.3.6'
    String vt = '/opt/ngenebio/app/vt/vt-0.5772-d00d603a/vt'

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
            log_progress 'Already finished to edit consecutive varinat using inhouse-script'
        else
            log_progress 'Start to edit consecutive varinat using inhouse-script'

            wrap python ${script} -i ${inputVcfFile} -o ${outputFile} -s ${outputName} -d ${outputDir} -b ${inputBamFile} -m ${consecutiveDb} -c ${freebayes} -r ${referenceFasta} -t ${vt} -f ${minAltFrac}

            wrap md5sum ${outputFile} > '${outputFile}.md5'
            log_progress 'Finish to edit consecutive varinat using inhouse-script'
        fi
    >>>

    output {
        String vcfFile = outputFile
    }
}