#================================================================================
#
# FILE: WhitelistFilter.wdl
#
# DESCRIPTION: NGeneBio DNA ONCOaccuPanel Whitelist Filter workflow
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
# VALIDATION : /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/womtool-50.jar validate /opt/ngenebio/pipeline/step_modules/SnvCalling/WhitelistProcess.wdl
#================================================================================

task makeTarget {

    # Input Parameter
    String stepName
    String sampleName
    String logFile
    String inputDisease
    String tmpDisease = sub(inputDisease , '_NOS$' , '')
    String inputTarget 
    String outputDir

    String outputTarget = outputDir + sampleName + '.whitelist.list.txt'
    
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
                
        ret=$(md5sum_check ${outputTarget})
        if [[ $ret = 'True' ]]; then
            log_progress 'Extract Diagnosis Target already finished!!!'
        else
            log_progress 'Extract Diagnosis Target start'

            if [[ ${inputDisease} == 'Cancer_Unknown_Primary_NOS' ]]; then
                cp ${inputTarget} ${outputTarget}

            elif [[ ${inputDisease} == 'CNS_Brain_Cancer_NOS' ]]; then
                grep -w 'Glioma\|Malignant_Brain_Tumor' ${inputTarget} > ${outputTarget}

            elif [[ ${inputDisease} == 'Lung_Cancer_NOS' ]]; then
                grep -w 'Small_Cell_Lung_Cancer\|Non_Small_Cell_Lung_Cancer' ${inputTarget} > ${outputTarget}

            elif [[ ${inputDisease} == 'Multiple_Myeloma_NOS' ]]; then
                grep -w 'Melanoma' ${inputTarget} > ${outputTarget}

            elif [[ ${inputDisease} == 'Soft_Tissue_Cancer_Sarcoma_NOS' ]]; then
                grep -w 'Gastrointestinal_Stromal_Tumor' ${inputTarget} > ${outputTarget}

            elif [[ ${inputDisease} == 'Uterine_Cancer_NOS' ]]; then
                grep -w 'Endometrial_Cancer' ${inputTarget} > ${outputTarget}

            else
                grep -w ${tmpDisease} ${inputTarget} > ${outputTarget}
            fi

            wrap md5sum ${outputTarget} > '${outputTarget}.md5'
            log_progress 'Extract Diagnosis Target finished'

        fi

    >>>

        output {
            String listFile = outputTarget
        }

}

task checkTarget {

    # Input Parameter
    String stepName
    String inputTarget
    String outputVcf
    String logFile

    String outputVariant = sub(outputVcf, '\\.final\\.vcf$','.variant.txt')

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

        isWhitelistTarget='true'

        if ! [[ -s ${inputTarget} ]]; then

            log_progress 'The disease does not have whitelist variants'
             > ${outputVariant}
            wrap md5sum ${outputVariant} > '${outputVariant}.md5'

            > ${outputVcf}
            
            wrap md5sum ${outputVcf} > '${outputVcf}.md5'

            isWhitelistTarget='false'

        fi
            echo $isWhitelistTarget

    >>>
        output{
            String finalTxt = read_string(stdout())
            String finalVariantFile = outputVariant
            String finalVcfFile = outputVcf
        }
}


task whitelistFilter {

    # Input Parameter
    String stepName
    String logFile
    String inputVcf #
    String inputTarget  ## output of MakeWhitelistTarget

    String outputVcf = sub(inputVcf, '\\.raw\\.vcf$','.filter.vcf') 
    String outputTxt = sub(inputVcf, '\\.raw\\.vcf$','.variant_tmp.txt') 

    # Used tools
    String script = '/opt/ngenebio/pipeline/step_modules/SnvCalling/onco_snvWhiFilter.py'

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

        ret1=$(md5sum_check ${outputVcf})
        ret2=$(vcf_check ${inputVcf})
        if [[ $ret1 = 'True' ]]; then
            log_progress 'Whitelist Filter already finished!!!'
        elif [[ $ret2 = 'False' ]]; then
            log_progress 'Whitelist Filter input vcf is empty!!!'

            wrap > ${outputVcf}
            wrap md5sum ${outputVcf} > '${outputVcf}.md5'
            log_progress 'Whitelist Filter finished'
        else
            log_progress 'Whitelist Filter start'

            wrap source /snv-pipeline/bin/activate

            wrap python ${script} -i ${inputVcf} -o ${outputVcf} -u ${outputTxt} -t ${inputTarget}

            wrap md5sum ${outputVcf} > '${outputVcf}.md5'
            log_progress 'Whitelist Filter finished'

            wrap deactivate
        fi
    >>>

    output {
        String vcfFile = outputVcf
        String txtFile = outputTxt
    }


}

task filterDepth {

    # Input Parameter
    String stepName
    String logFile
    String inputVcf ## whiltelist_format_vcf
    String vcfFormat 
    Int cutoff

    String outputVcf  = sub(inputVcf, '\\.reform\\.vcf','.depth.vcf')  

    # Used tools
    String script = '/opt/ngenebio/pipeline/step_modules/SnvCalling/onco_snvFilter.py'
    String snvRecover = '/opt/ngenebio//pipeline/assay_reference/ONCO/COMMON/recoverVariant.json'
    
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
        ret1=$(md5sum_check ${outputVcf})
        ret2=$(vcf_check ${inputVcf})

        if [[ $ret1 = 'True' ]]; then
            log_progress 'Filter vcf depth already finished!!!'
        elif [[ $ret2 = 'False' ]]; then
            log_progress 'Filter vcf depth input vcf is empty!!!'

            wrap > ${outputVcf}
            wrap md5sum ${outputVcf} > '${outputVcf}.md5'
            log_progress 'Filter vcf depth finished'
        else
            log_progress 'Filter vcf depth start'

            wrap source /snv-pipeline/bin/activate
            wrap python ${script} --input ${inputVcf} --output ${outputVcf} --format ${vcfFormat} \
                                  --cutoff ${cutoff} --recover ${snvRecover}

            wrap md5sum ${outputVcf} > '${outputVcf}.md5'
            log_progress 'Filter vcf depth finished'

            wrap deactivate
        fi

    >>>

    output {
        String vcfFile = outputVcf
    }
}

task whitelistInfo {

    # Input Parameter
    String stepName
    String logFile
    String inputVcf ## whitelist Final VCF
    String inputTxt  ## whitelistFilter output Txt

    String outputTxt = sub(inputVcf, '\\.final\\.vcf$','.variant.txt')  
    
    # Used tools
    String script = '/opt/ngenebio/pipeline/step_modules/SnvCalling/onco_snvWhiInfo.py'

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

        ret1=$(md5sum_check ${outputTxt})
        ret2=$(vcf_check ${inputVcf})
        if [[ $ret1 = 'True' ]]; then
            log_progress 'Set whitelist variant info already finished!!!'
        elif [[ $ret2 = 'False' ]]; then
            log_progress 'Set whitelist variant info input vcf is empty!!!'

            wrap > ${outputTxt}
            wrap md5sum ${outputTxt} > '${outputTxt}.md5'
            log_progress 'Set whitelist variant info finished'
        else
            wrap source /snv-pipeline/bin/activate
            log_progress 'Set whitelist variant info start'
            wrap python ${script} -i ${inputVcf} -j ${inputTxt} -o ${outputTxt}

            wrap md5sum ${outputTxt} > '${outputTxt}.md5'
            log_progress 'Set whitelist variant info finished'

            wrap deactivate
        fi
    >>>
    output {
        String txtFile = outputTxt
    }

}

# task vep {

#     String inputVcf  
#     String inputFasta
#     String stepName = 'vep'
#     String logFile 
#     String outputVcf   
#     String tmpVep = sub(outputVcf, '\\.vcf$','_tmp.vep.vcf')
#     String vep = '/opt/ngenebio/app/vep/variant_effect_predictor.pl'
#     String vepPath = '/opt/ngenebio/app/vep/'
#     String vepConverter = '/opt/ngenebio/utils/vepconverter/snv_vep_converter.py'
#     String enstamc = '/opt/ngenebio/pipeline/assay_reference/ONCO/V1/isoform_overrides_amc'

#     command <<<
#         export inputStep=${stepName}
#         export logFile=${logFile}
#         source /opt/ngenebio/pipeline/script/validate.sh

#         ret1=$(md5sum_check ${outputVcf})
#         ret2=$(vcf_check ${inputVcf})
#         if [[ $ret1 = 'True' ]]; then
#             log_progress 'VEP annotation already finished!!!'
#         elif [[ $ret2 = 'False' ]]; then
#             log_progress 'VEP annotation input vcf is empty!!!'

#             wrap > ${outputVcf}
#             wrap md5sum ${outputVcf} > '${outputVcf}.md5'
#             log_progress 'VEP annotation finished'
#         else
#             log_progress 'VEP annotation start'

#             wrap source /annotation-pipeline/bin/activate
#             wrap perl ${vep} --species homo_sapiens --assembly GRCh37 --no_progress --no_stats --buffer_size 5000 --sift b    \
#                                 --ccds --uniprot --hgvs --symbol --numbers --domains --gene_phenotype --canonical --protein      \
#                                 --biotype --uniprot --tsl --variant_class --shift_hgvs 1 --check_existing --total_length         \
#                                 --allele_number --no_escape --xref_refseq --failed 1 --vcf --flag_pick_allele --pick_order       \
#                                 canonical,tsl,biotype,rank,ccds,length --dir ${vepPath} --fasta ${inputFasta} --format vcf   \
#                                 --offline --cache_version 92 --force_overwrite --input_file ${inputVcf} --output_file ${tmpVep}

#             wrap python ${vepConverter} --input ${tmpVep} --output ${outputVcf} --custom_enst ${enstamc}

#             wrap md5sum ${outputVcf} > '${outputVcf}.md5'
#             log_progress 'VEP annotation finished'

#             wrap deactivate
#         fi

#     >>>
#     output {
#         String finalVcfFile = outputVcf
#     }

# }