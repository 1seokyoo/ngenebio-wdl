#================================================================================
#
# FILE: VcfMerge.wdl
#
# DESCRIPTION: NGeneBio DNA ONCOaccuPanel VCF Merge workflow
# Note:
# REQUIREMENTS:
# BUGS:
# NOTES:
# AUTHOR: Juyeon hong hong.juyeon@ngenebio.com
# COMPANY: ngenebio
# VERSION: 0.1
# CREATED: 2022.05.17
# REVISION:
# RUN COMMAND:
# VALIDATION : /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/womtool-50.jar validate /opt/ngenebio/pipeline/step_modules/SnvCalling/VcfMerge.wdl
#================================================================================

task vcfMerge{

    # Input Parameter
    String stepName 
    String inputAllGeneVcf
    String inputHrdVcf
    String inputWhitelistVcf
    String inputMutectVcf
    String inputSnvFormat  
    String sampleName
    String outputDir
    String logFile

    String outputVcf
    String tmpHdrFile = outputDir + sampleName + '.total.raw.hdr'
    String tmpTotalRawVcf = outputDir + sampleName + '.total.tmp.raw.vcf'
    String tmpReformVcf = outputDir + sampleName + '.total.tmp.reform.vcf'

    # Used tools
    String script = '/opt/ngenebio/pipeline/step_modules/SnvCalling/onco_snvMerge.py'
    String vt = '/opt/ngenebio/app/vt/vt-0.5772-d00d603a/vt'
    
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

        ret=$(md5sum_check ${outputVcf})
        if [[ $ret = 'True' ]]; then
            log_progress 'Merge vcf already finished!!!'
        else
            log_progress 'Merge vcf start'

            wrap source /snv-pipeline/bin/activate

            wrap sed -e 's%SAMPLE_NAME%${sampleName}%' ${inputSnvFormat} > ${tmpHdrFile}

            wrap python ${script} -n ${inputAllGeneVcf} -w ${inputWhitelistVcf} -r ${inputHrdVcf} \
                                  -m ${inputMutectVcf} -o ${tmpTotalRawVcf} -t ${tmpHdrFile}

            wrap md5sum ${tmpTotalRawVcf} > '${tmpTotalRawVcf}.md5'

            grep -v '<DUP>\|<DEL>\|<INV>' ${tmpTotalRawVcf} > ${tmpReformVcf}

            wrap ${vt} sort -o ${outputVcf} ${tmpReformVcf}

            wrap md5sum ${outputVcf} > '${outputVcf}.md5'
            log_progress 'Merge vcf finished'

            wrap deactivate
        fi
    >>>

    output {
        String vcfFile = outputVcf
    }
}