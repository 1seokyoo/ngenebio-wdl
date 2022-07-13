#================================================================================
#
# FILE: RunMutect2.wdl
#
# DESCRIPTION: NGeneBio DNA ONCOaccuPanel GATK 4.2~~ Mutect2 workflow
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
# VALIDATION : /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/womtool-50.jar validate /opt/ngenebio/pipeline/step_modules/SnvCalling/RunMutect2.wdl
#================================================================================

task cleanSam {

    # Input Parameter
    String stepName
    String inputBam #fianl_bam_out
    String outputDir 
    String sampleName
    String logFile
    String tmpDir
    Int cpu
    Int ram

    String outputBam = outputDir + sampleName + '.mutect.clean.bam'
    String tmpBam = sub(outputBam, '\\.bam$','_tmp.bam')

    # Used tools
    String gatk = '/opt/ngenebio/app/GATK/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar'
    String java = '/opt/ngenebio/app/java/jre1.8.0_65/bin/java'
    String samtools = '/opt/ngenebio/app/samtools/samtools-0.1.19/samtools'

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

        ret=$(md5sum_check ${outputBam})
        if [[ $ret = 'True' ]]; then
            log_progress 'Cleans final BAM already finished!!!'
        else
            log_progress 'Cleans final BAM start'
            
            wrap source /snv-pipeline/bin/activate
            wrap ${samtools} view -b -F 4 -q 1 -b -h -o ${tmpBam} ${inputBam}

            wrap ${java} -Xmx${ram}g -XX:ParallelGCThreads=${cpu} -Djava.io.tmpdir=${tmpDir} \
                         -jar ${gatk} CleanSam \
                                      --INPUT ${tmpBam} \
                                      --OUTPUT ${outputBam} \
                                      --CREATE_INDEX true

            wrap md5sum ${outputBam} > '${outputBam}.md5'
            log_progress 'Cleans final BAM finished'

            wrap deactivate
        fi
    >>>

    output {
        String bamFile = outputBam
    }
}

task mutect2 {
    
    # Input Parameter
    String stepName
    String inputBam 
    String inputTarget 
    String inputFasta
    String sampleName
    String outputDir
    String logFile
    String tmpDir
    Int cpu
    Int ram

    String outputVcf = outputDir + sampleName + '.mutect.raw.vcf.gz'
    String outputBam = outputDir + sampleName + '.mutect.bam'
    String outputRegionTxt = outputDir + sampleName + '.assembly_region.txt'
    String outputF1R1 = outputDir + sampleName + '.f1r2.tar.gz'

    # Used tools
    String gatk = '/opt/ngenebio/app/GATK/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar'
    String java = '/opt/ngenebio/app/java/jre1.8.0_65/bin/java'
    String samtools = '/opt/ngenebio/app/samtools/samtools-0.1.19/samtools'

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
        ret2=$(md5sum_check ${outputF1R1})
        if [[ $ret1 = 'True' && $ret2 = 'True' ]]; then
            log_progress 'Mutect2 variant calling already finished!!!'
        else
            log_progress 'Mutect2 variant calling start'

            wrap source /snv-pipeline/bin/activate
            wrap ${java} -Xmx${ram}g -XX:ParallelGCThreads=${cpu} -Djava.io.tmpdir=${tmpDir} \
                         -jar ${gatk} Mutect2 --input ${inputBam} --output ${outputVcf} \
                                              --reference ${inputFasta} --tumor-sample ${sampleName} \
                                              --bam-output ${outputBam} --assembly-region-out ${outputRegionTxt} \
                                              --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
                                              --intervals ${inputTarget} --native-pair-hmm-threads ${cpu} \
                                              --max-mnp-distance 3 --f1r2-tar-gz ${outputF1R1}

            ### Mutect2 BAM indexing
            log_progress 'Mutect2 BAM indexing start'
            wrap ${samtools} index ${outputBam}

            wrap md5sum ${outputVcf} > '${outputVcf}.md5'
            wrap md5sum ${outputF1R1} > '${outputF1R1}.md5'

            log_progress 'Mutect2 variant calling finished'

            wrap deactivate
        fi
    >>>

    output {
            String bamFile = outputBam
            String vcfFile = outputVcf
            String f1R2File = outputF1R1
    }    
}


task mutect2PostProcessing{

    # Input Parameter
    String stepName
    String inputBam
    String inputVcf
    String inputTarget
    String inputFasta
    String inputF1R2
    String sampleName
    String outputDir
    String tmpDir
    String logFile
    Int cpu
    Int ram

    String outputVcf = outputDir + sampleName + '.mutect.post.vcf'
    String tmpPileupSummary = sub(outputVcf, '\\.vcf$','.pileupsummary.table')
    String tmpCalculateContam = sub(outputVcf, '\\.vcf$','.calculatecontam.table')
    String tmpOrientationModel = outputDir + sampleName  + '.read.orientation.model.tar.gz'
    
    # Used tools
    String gatk = '/opt/ngenebio/app/GATK/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar'
    String java = '/opt/ngenebio/app/java/jre1.8.0_65/bin/java'
    String ExAC = '/opt/ngenebio/dependencies/ExAC/ExAC.hg19.vcf.gz'

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
            log_progress 'Mutect2 post-processing already finished!!!'
        else

            log_progress 'MuTect2 post-processing start'
            wrap source /snv-pipeline/bin/activate

            ## Run LearnReadOrientationModel
            log_progress 'LearnReadOrientationModel start'
            wrap ${java} -Xmx${ram}g -XX:ParallelGCThreads=${cpu} -Djava.io.tmpdir=${tmpDir}   \
                -jar ${gatk} LearnReadOrientationModel -I ${inputF1R2} -O ${tmpOrientationModel}
            log_progress 'LearnReadOrientationModel done'

            ### Run Pileup Summary ###
            log_progress 'GetPileupSummaries start'
            wrap ${java} -Xmx${ram}g -XX:ParallelGCThreads=${cpu} -Djava.io.tmpdir=${tmpDir}   \
                -jar ${gatk} GetPileupSummaries --input ${inputBam}                            \
                                                --output ${tmpPileupSummary}                   \
                                                --variant ${ExAC} --intervals ${inputTarget}
            log_progress 'GetPileupSummaries done'

            ### Run Calculate Contamination ###
            log_progress 'CalculateContamination start'
            wrap ${java} -Xmx${ram}g -XX:ParallelGCThreads=${cpu} -Djava.io.tmpdir=${tmpDir}   \
                -jar ${gatk} CalculateContamination --input ${tmpPileupSummary}                \
                                                    --output ${tmpCalculateContam}
            log_progress 'CalculateContamination done'

            ### Run Filter Mutect Calls ###
            log_progress 'FilterMutectCalls start'
            wrap ${java} -Xmx${ram}g -XX:ParallelGCThreads=${cpu} -Djava.io.tmpdir=${tmpDir}    \
                -jar ${gatk} FilterMutectCalls  --variant ${inputVcf}                           \
                                                --output ${outputVcf}                           \
                                                --contamination-table ${tmpCalculateContam}     \
                                                --ob-priors ${tmpOrientationModel}              \
                                                --reference ${inputFasta}
            log_progress 'FilterMutectCalls done'

            wrap md5sum ${outputVcf} > '${outputVcf}.md5'
            log_progress 'MuTect2 post-processing done'

            wrap deactivate
        fi
    >>>

    output{
        String vcfFile = outputVcf
    }
}