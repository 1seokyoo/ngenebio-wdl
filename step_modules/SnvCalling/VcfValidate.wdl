#================================================================================
#
# FILE: MakeWhitelistTarget.wdl
#
# DESCRIPTION: NGeneBio DNA ONCOaccuPanel Make Whitelist Target workflow
# Note:
# REQUIREMENTS:
# BUGS:
# NOTES:
# AUTHOR: Juyeon hong hong.juyeon@ngenebio.com
# COMPANY: ngenebio
# VERSION: 0.1
# CREATED: 2022.04.24
# REVISION:
# RUN COMMAND: 
# VALIDATION : /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/womtool-50.jar validate /opt/ngenebio/pipeline/step_modules/SnvCalling/VcfValidate.wdl
#================================================================================

task vcfValidate{

    # Input Parameter
    String stepName
    String inputVcf  
    String inputFasta
    String sampleName
    String outputDir
    String mode
    String tmpDir
    String logFile
    Int cpu
    Int ram

    Boolean? active_override # true or false
    Boolean active_value = select_first([active_override, false])
    String active_type = if active_value then 'True' else 'False'

    String outputVcf  = outputDir + sampleName + '.' + mode + '.reform.vcf'
    String outputSortedVcf = sub(outputVcf, '\\.vcf$','_tmp.sorted.vcf')
    String outputConfig = sub(outputVcf, '\\.vcf$','_tmp.config.vcf')
    String outputDecomposeVcf = sub(outputVcf, '\\.vcf$','_tmp.decompose.vcf')
    String outputNormalizeVcf = sub(outputVcf, '\\.vcf$','_tmp.normalize.vcf')

    # Used tools
    String java = '/opt/ngenebio/app/java/jre1.8.0_65/bin/java'
    String vt = '/opt/ngenebio/app/vt/vt-0.5772-d00d603a/vt'
    String gatk = '/opt/ngenebio/app/GATK/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar'

    command <<<

        export inputStep=${stepName}
        export logFile=${logFile}
        source /opt/ngenebio/pipeline/script/validate.sh

        ret1=$(md5sum_check ${outputVcf})
        ret2=$(vcf_check ${inputVcf})
        if [[ $ret1 = 'True' ]]; then
            log_progress 'VCF Validation already finished!!!'
        elif [[ $ret2 = 'False' ]]; then
            log_progress 'VCF Validation input vcf is empty!!!'

            wrap > ${outputVcf}
            wrap md5sum ${outputVcf} > '${outputVcf}.md5'
            log_progress 'VCF Validation finished'
        else
            log_progress 'VCF Validation start'
            log_progress ${active_value}
            log_progress ${active_type}

            if [[ ${active_type} = 'True' ]]; then
                wrap source /snv-pipeline/bin/activate
            fi 

            wrap ${java} -Xmx${ram}g -XX:ParallelGCThreads=${cpu} -Djava.io.tmpdir=${tmpDir}    \
                            -jar ${gatk} SelectVariants -V ${inputVcf} -R ${inputFasta} -O ${outputConfig}

            wrap ${vt} sort -o ${outputSortedVcf} ${outputConfig}

            wrap ${vt} decompose -s -o ${outputDecomposeVcf} ${outputSortedVcf}

            wrap ${vt} normalize -r ${inputFasta} -o ${outputNormalizeVcf} ${outputDecomposeVcf}

            wrap cp ${outputNormalizeVcf} ${outputVcf}

            wrap md5sum ${outputVcf} > '${outputVcf}.md5'
            log_progress 'VCF Validation finished'

            if [[ ${active_type} = 'True' ]]; then
                deactivate
            fi 
        fi

    >>>

    output {
        String formatVcfFile = outputVcf
    }

}

