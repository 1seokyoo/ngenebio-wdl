#================================================================================
#
# FILE: ONCOBwaRun.wdl
#
# DESCRIPTION: NGeneBio DNA fastq alignment pipeline
# Note:
# REQUIREMENTS:
# BUGS:
# NOTES:
# AUTHOR: Juyeon Hong hong.juyeon@ngenebio.com
# COMPANY: ngenebio
# VERSION: 0.1
# CREATED: 2022.04.25
# REVISION:
# RUN_COMMAND: /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/cromwell-50.jar run /opt/ngenebio/pipeline/step_modules/DnaAlignment/BwaV059.wdl -i /opt/ngenebio/pipeline/step_modules/DnaAlignment/testBwaV059Run.json
# VALIDATION: /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/womtool-50.jar validate /opt/ngenebio/pipeline/step_modules/DnaAlignment/BwaV059.wdl
#================================================================================

task bwaAln{

    # Input Parameter
    String stepName
    String readStrand
    String logFile
    String inputFasta
    String inputFastq
    String outputDir
    String bwa = '/opt/ngenebio/app/BWA/bwa-0.5.9/bwa'
    Int cpu
    
    # Set Parameter
    String outputSai = outputDir + sub(basename(inputFastq), '\\.fastq.gz$', '.raw.sai')  
    
    command <<<
        export inputStep=${stepName}
        export logFile=${logFile}
        source /opt/ngenebio/pipeline/script/validate.sh

        ret=$(md5sum_check ${outputSai})
        if [[ $ret = 'True' ]]; then
                log_progress 'Fastq indexing already finished!!!'
            else
                log_progress 'Fastq indexing start'
                (wrap ${bwa} aln -t ${cpu} -q 5 -l 32 -k 2 -o 1 -f ${outputSai} ${inputFasta} ${inputFastq}) 
                wrap md5sum ${outputSai} > ${outputSai}.md5

                log_progress 'Fastq indexing finished'
            fi

    >>>

    output {
        String saiFile = outputSai
    }
}

task bwaSampe{

    # Input Parameter
    String stepName
    String logFile
    String inputFastqR1
    String inputFastqR2
    String inputSaiR1
    String inputSaiR2
    String inputFasta
    String outputDir
    String sampleName
    String tmpDir 
    String bwa = '/opt/ngenebio/app/BWA/bwa-0.5.9/bwa'
    String samtools = '/opt/ngenebio/app/samtools/samtools-0.1.19/samtools'
    String java = '/opt/ngenebio/app/java/jre1.8.0_65/bin/java'
    String picard = '/opt/ngenebio/app/picard/picard-tools-2.5.0/picard.jar'
    Int cpu
    Int ram


    # Set Parameter
    String outputFile = outputDir + sampleName + '.raw.bam' ## outputDir 'rundir/data/basecall/alignment/'
    String tmpOutputBam =  outputDir + sampleName + '.nosorted.bam'  ## outputDir 'rundir/data/basecall/alignment/'
    String tmpOutputSam = outputDir + sampleName + '.nosorted.sam' ## outputDir 'rundir/data/basecall/alignment/'
     
    command <<<
        export inputStep=${stepName}
        export logFile=${logFile}
        source /opt/ngenebio/pipeline/script/validate.sh

        ret=$(md5sum_check ${outputFile})
        if [[ $ret = 'True' ]]; then
            log_progress 'BWA alignment & Picard sort already finished!!!'
        else
            log_progress 'BWA alignment & Picard sort start'
            wrap ${bwa} sampe -f ${tmpOutputSam} -P ${inputFasta} ${inputSaiR1} ${inputSaiR2} ${inputFastqR1} ${inputFastqR2} \
                -r "'@RG\\tID:${sampleName}\\tPL:illumina\\tPU:ex\\tLB:${sampleName}\\tSM:${sampleName}'"

            wrap ${samtools} view -b -h -S -o ${tmpOutputBam} ${tmpOutputSam}
            wrap ${samtools} index ${tmpOutputBam}
            wrap rm -f ${tmpOutputSam}
            wrap rm -f ${tmpOutputSam}.bai
            log_progress 'BWA alignment finished'

            log_progress 'Picard sort start'
            wrap ${java} -XX:ParallelGCThreads=${cpu} -Xms${ram}g -Xmx${ram}g                        \
                -Dsamjdk.compression_level=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10                        \
                -jar ${picard} SortSam INPUT=${tmpOutputBam} OUTPUT=${outputFile} TMP_DIR=${tmpDir}    \
                SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true

            wrap ${samtools} index ${outputFile}
            wrap md5sum ${outputFile} > '${outputFile}.md5'
            wrap rm -f ${tmpOutputBam}
            log_progress 'BWA alignment & Picard sort finished'
        fi

    >>>

    output {
        String bamFile = outputFile
    }
}
