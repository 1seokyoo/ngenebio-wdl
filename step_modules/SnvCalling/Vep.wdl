#================================================================================
#
# FILE: Vep.wdl
#
# DESCRIPTION: NGeneBio DNA ONCOaccuPanel VEP workflow
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
# VALIDATION : /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/womtool-50.jar validate /opt/ngenebio/pipeline/step_modules/SnvCalling/Vep.wdl
#================================================================================


task run {

    # Input Parameter
    String stepName
    String inputVcf
    String inputFasta
    String sampleName
    String outputDir
    Int ensemblVersion
    String logFile

    String outputVcf = outputDir + sampleName + 'total.vep.vcf'
    String tmpVepVcf = outputDir + sampleName + 'total.vep.tmp.vcf'
    
    # Used tools
    String vep = '/opt/ngenebio/app/vep/variant_effect_predictor.pl'
    String vepPath = '/opt/ngenebio/app/vep/'    

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
            log_progress 'VEP annotation already finished!!!'
        else
            log_progress 'VEP annotation start'
            wrap source /annotation-pipeline/bin/activate

            wrap perl ${vep} --species homo_sapiens --assembly GRCh37 --no_progress --no_stats --buffer_size 5000 \
                             --sift b --ccds --uniprot --hgvs --symbol --numbers --domains --gene_phenotype --canonical \
                             --protein --biotype --uniprot --tsl --variant_class --shift_hgvs 1 --check_existing \
                             --total_length --allele_number --no_escape --xref_refseq --failed 1 --vcf --flag_pick_allele \
                             --pick_order canonical,tsl,biotype,rank,ccds,length --dir ${vepPath} --fasta ${inputFasta} \
                             --format vcf --offline --cache_version ${ensemblVersion} --force_overwrite \
                             --input_file ${inputVcf} --output_file ${outputVcf}

            wrap md5sum ${outputVcf} > '${outputVcf}.md5'
            log_progress 'VEP annotation finished'
            wrap deactivate
        fi
    >>>

    output {
        String vcfFile = outputVcf
    }
}


task filter {

    # Input Parameter
    String stepName
    String inputVcf
    String sampleName
    String outputDir
    String logFile

    String outputVcf = outputDir + sampleName + 'total.filter.vcf'
    String tmpVepVcf = outputDir + sampleName + 'total.filter.tmp.vcf'
    
    # Used tools
    String script = '/opt/ngenebio/pipeline/step_modules/SnvCalling/onco_snvVepConverter.py'
    String enstamc = '/opt/ngenebio/pipeline/assay_reference/ONCO/COMMON/isoform_overrides_amc'

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
            log_progress 'SNV/INDEL VEP filter already finished!!!'
        else
            log_progress 'SNV/INDEL VEP filter start'
            wrap source /annotation-pipeline/bin/activate

            wrap python ${script} --input ${inputVcf} --output ${tmpVepVcf} --custom_enst ${enstamc}

            awk '$7 != "REJECT" && $7 != "Mutect2"' ${tmpVepVcf} > ${outputVcf}

            wrap md5sum ${outputVcf} > '${outputVcf}.md5'
            log_progress 'SNV/INDEL VEP filter finished'
            wrap deactivate
        fi

    >>>

    output {
        String vcfFile = outputVcf 
    }
}



