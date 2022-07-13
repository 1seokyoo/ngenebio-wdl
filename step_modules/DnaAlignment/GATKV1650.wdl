#================================================================================
#
# FILE: GATK.wdl
#
# DESCRIPTION: NGeneBio DNA GATK pipeline
# Note:
# REQUIREMENTS:
# BUGS:
# NOTES:
# AUTHOR: Juyeon Hong hong.juyeon@ngenebio.com
# COMPANY: ngenebio
# VERSION: 0.1
# CREATED: 2022.04.25
# REVISION:
# RUN_COMMAND:
# VALIDATION: /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/womtool-50.jar validate /opt/ngenebio/pipeline/step_modules/DnaAlignment/GATKV1650.wdl
#================================================================================


task indelRealign{

    # Input Parameter
    String stepName
    String inputBam
    String inputFasta
    String inputIntervalBaits
    String inputMills
    String inputG1000Phase1
    String logFile
    String tmpDir
    Int cpu
    Int ram

    String outputBam = sub(inputBam, '\\.dp.bam$', '.realign.bam')
    
    # Used tools
    String java = '/opt/ngenebio/app/java/jre1.7.0_60/bin/java'
    String gatk = '/opt/ngenebio/app/GATK/gatk-1.6.5.0/GenomeAnalysisTK.jar'
    
    command <<<
        export inputStep=${stepName}
        export logFile=${logFile}
        source /opt/ngenebio/pipeline/script/validate.sh

        ret=$(md5sum_check ${outputBam})
        if [[ $ret = 'True' ]]; then
            log_progress 'Indel Realign already finished!!!'
        else
            log_progress 'Indel Realign start'
            wrap ${java} -XX:ParallelGCThreads=${cpu}  -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 \
                         -Xms${ram}g -Xmx${ram}g -Djava.io.tmpdir=${tmpDir} \
                         -jar ${gatk} -T IndelRealigner \
                                      -I ${inputBam} \
                                      -o ${outputBam} \
                                      -R ${inputFasta} \
                                      -targetIntervals ${inputIntervalBaits} \
                                      -known ${inputMills} \
                                      -known ${inputG1000Phase1} \
                                      -model KNOWNS_ONLY \
                                      -LOD 0.4 \
                                      -compress 1 \
                                      -maxInMemory 1000000                                   

            grep -v 'WARN' ${logFile}
            wrap md5sum ${outputBam} > '${outputBam}.md5'
            log_progress 'Indel Realign finished'
        fi

    >>>

    output {
        String realignBamFile = outputBam
    }

}


task countCovariates{

    # Input Parameter
    String stepName
    String inputBam
    String inputFasta
    String dbsnp 
    String logFile
    String tmpDir
    Int cpu
    Int ram

    String outputCsv = sub(inputBam, '\\.realign.bam$', '.covariate.csv')
    
    # Used tools
    String java = '/opt/ngenebio/app/java/jre1.7.0_60/bin/java'
    String gatk = '/opt/ngenebio/app/GATK/gatk-1.6.5.0/GenomeAnalysisTK.jar'
    
    command <<<
        export inputStep=${stepName}
        export logFile=${logFile}
        source /opt/ngenebio/pipeline/script/validate.sh

        ret=$(md5sum_check ${outputCsv})
        if [[ $ret = 'True' ]]; then
            log_progress 'Count Covariates already finished!!!'
        else
            log_progress 'Count Covariates start'
            wrap ${java} -XX:ParallelGCThreads=${cpu}  -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 \
                         -Xms${ram}g -Xmx${ram}g -Djava.io.tmpdir=${tmpDir} \
                         -jar ${gatk} -T CountCovariates \
                                      -I ${inputBam} \
                                      -R ${inputFasta} \
                                      -recalFile ${outputCsv} \
                                      -knownSites ${dbsnp} \
                                      -l INFO \
                                      -nt 1 \
                                      -cov ReadGroupCovariate \
                                      -cov CycleCovariate \
                                      -cov DinucCovariate \
                                      -cov QualityScoreCovariate \
                                      -OQ

            grep -v 'WARN' ${logFile}
            wrap md5sum ${outputCsv} > '${outputCsv}.md5'
            log_progress 'Count Covariates finished'
        fi

    >>>

    output {
        String covariateCsvFile = outputCsv
    }

}

task baseRecalibration{

    # Input Parameter
    String stepName
    String inputBam
    String inputCsv
    String inputFasta
    String logFile
    String tmpDir
    Int cpu
    Int ram

    String outputBam = sub(inputBam, '\\.realign.bam$', '.recal.bam')

    # Used tools
    String java = '/opt/ngenebio/app/java/jre1.7.0_60/bin/java'
    String gatk = '/opt/ngenebio/app/GATK/gatk-1.6.5.0/GenomeAnalysisTK.jar'
    String samtools = '/opt/ngenebio/app/samtools/samtools-0.1.19/samtools'
    
    command <<<
        export inputStep=${stepName}
        export logFile=${logFile}
        source /opt/ngenebio/pipeline/script/validate.sh

        ret=$(md5sum_check ${outputBam})
        if [[ $ret = 'True' ]]; then
            log_progress 'Base Recalibrator already finished!!!'
        else
            log_progress 'Base Recalibrator start'
            wrap ${java} -XX:ParallelGCThreads=${cpu}  -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 \
                         -Xms${ram}g -Xmx${ram}g -Djava.io.tmpdir=${tmpDir} \
                         -jar ${gatk} -T TableRecalibration \
                                      -I ${inputBam} \
                                      -R ${inputFasta} \
                                      --out ${outputBam} \
                                      -recalFile ${inputCsv} \
                                      -OQ                           

            grep -v 'WARN' ${logFile}
            wrap ${samtools} index ${outputBam}
            wrap md5sum ${outputBam} > '${outputBam}.md5'
            log_progress 'Base Recalibrator finished'
        fi
    >>>

    output {
        String recalBamFile = outputBam
    }

}