#================================================================================
#
# FILE: GATK.wdl
#
# DESCRIPTION: NGeneBio DNA GATK v4.2.6.1 pipeline
# Note:
# REQUIREMENTS:
# BUGS:
# NOTES:
# AUTHOR: Wonseok Yoo wonseok.yoo@ngenebio.com
# COMPANY: ngenebio
# VERSION: 0.1
# CREATED: 2022.05.28
# REVISION:
# RUN_COMMAND: 
# VALIDATION: 
#================================================================================

task addOrReplaceReadGroups{

    # Input Parameter
    String stepName
    String logFile
    String sampleName
    String inputFile
    String outputDir
    String outputName
    String tmpDir
    Int cpu
    Int memory

    String outputFile = outputDir + outputName + '.sort.bam'
    String rgpl = "Illumina"
    String rgpu = "NGeneBioBRCA1n2AMP"

    # Used tools
    String java = '/opt/ngenebio/app/java/jre1.8.0_65/bin/java'
    String gatk = '/opt/ngenebio/app/GATK/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar'
    String samtools = '/opt/ngenebio/app/samtools/samtools-1.3.1/samtools'

    command <<<

        export inputStep=${stepName}
        export logFile=${logFile}
        source /opt/ngenebio/pipeline/script/validate.sh

        ret=$(md5sum_check ${outputFile})
        if [[ $ret = 'True' ]]; then
            log_progress 'Already finished to base recalibrate using GATK v4.2.6.1'
        else
            log_progress 'Start to base recalibrate using GATK v4.2.6.1'

            wrap ${java} -Xmx${memory}g -XX:ParallelGCThreads=${cpu} -Djava.io.tmpdir=${tmpDir} \
                         -jar ${gatk} AddOrReplaceReadGroups \
                                        --INPUT ${inputFile} \
                                        --OUTPUT ${outputFile} \
                                        --RGLB ${sampleName} \
                                        --RGPL ${rgpl} \
                                        --RGPU ${rgpu} \
                                        --RGSM ${sampleName} \
                                        --RGID ${sampleName} \
                                        --SORT_ORDER coordinate \
                                        --VALIDATION_STRINGENCY LENIENT \
                                        2>&1
            wrap ${samtools} index ${outputFile}

            wrap md5sum ${outputFile} > '${outputFile}.md5'
            log_progress 'Finish to base recalibrate using GATK v4.2.6.1'
        fi
    >>>

    output {
        String bamFile = outputFile
    }
}

task baseRecalibration{

    # Input Parameter
    String stepName
    String logFile
    String sampleName
    String inputFile
    String outputDir
    String tmpDir
    String referenceFasta
    String dbsnp
    String mills
    String g1000Phase1
    Int cpu
    Int memory

    String outputFile = outputDir + sampleName + '.recal.grp'

    # Used tools
    String java = '/opt/ngenebio/app/java/jre1.8.0_65/bin/java'
    String gatk = '/opt/ngenebio/app/GATK/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar'

    command <<<

        export inputStep=${stepName}
        export logFile=${logFile}
        source /opt/ngenebio/pipeline/script/validate.sh

        ret=$(md5sum_check ${outputFile})
        if [[ $ret = 'True' ]]; then
            log_progress 'Already finished to base recalibrate using GATK v4.2.6.1'
        else
            log_progress 'Start to base recalibrate using GATK v4.2.6.1'

            wrap ${java} -Xmx${memory}g -XX:ParallelGCThreads=${cpu} -Djava.io.tmpdir=${tmpDir} \
                         -jar ${gatk} BaseRecalibrator \
                                        --input ${inputFile} \
                                        --output ${outputFile} \
                                        --reference ${referenceFasta} \
                                        --known-sites ${dbsnp} \
                                        --known-sites ${mills} \
                                        --known-sites ${g1000Phase1}

            wrap md5sum ${outputFile} > '${outputFile}.md5'
            log_progress 'Finish to base recalibrate using GATK v4.2.6.1'
        fi
    >>>

    output {
        String recalFile = outputFile
    }
}

task applyBqsr{

    # Input Parameter
    String stepName
    String logFile
    String sampleName
    String inputBamFile
    String inputRecal
    String outputDir
    String tmpDir
    String referenceFasta
    Int cpu
    Int memory

    String outputFile = outputDir + sampleName + '.bqsr.bam'

    # Used tools
    String java = '/opt/ngenebio/app/java/jre1.8.0_65/bin/java'
    String gatk = '/opt/ngenebio/app/GATK/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar'
    String samtools = '/opt/ngenebio/app/samtools/samtools-1.3.1/samtools'

    command <<<

        export inputStep=${stepName}
        export logFile=${logFile}
        source /opt/ngenebio/pipeline/script/validate.sh

        ret=$(md5sum_check ${outputFile})
        if [[ $ret = 'True' ]]; then
            log_progress 'Already finished to apply BQSR using GATK v4.2.6.1'
        else
            log_progress 'Start to apply BQSR using GATK v4.2.6.1'

            wrap ${java} -Xmx${memory}g -XX:ParallelGCThreads=${cpu} -Djava.io.tmpdir=${tmpDir} \
                         -jar ${gatk} ApplyBQSR \
                                        --input ${inputBamFile} \
                                        --bqsr-recal-file ${inputRecal} \
                                        --output ${outputFile} \
                                        --reference ${referenceFasta}

            wrap ${samtools} index ${outputFile}

            wrap md5sum ${outputFile} > '${outputFile}.md5'
            log_progress 'Finish to apply BQSR using GATK v4.2.6.1'
        fi
    >>>

    output {
        String bamFile = outputFile
    }
}