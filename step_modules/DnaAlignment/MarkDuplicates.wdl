#================================================================================
#
# FILE: MarkDuplicates.wdl
#
# DESCRIPTION: NGeneBio DNA MarkDuplicates pipeline
# Note:
# REQUIREMENTS:
# BUGS:
# NOTES:
# AUTHOR: Juyeon Hong hong.juyeon@ngenebio.com
# COMPANY: ngenebio
# VERSION: 0.1
# CREATED: 2022.04.25
# REVISION:
# RUN_COMMAND: /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/cromwell-50.jar run /opt/ngenebio/pipeline/step_modules/DnaAlignment/MarkDuplicates.wdl -i /opt/ngenebio/pipeline/step_modules/DnaAlignment/testMarkDuplicates.json
# VALIDATION: /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/womtool-50.jar validate /opt/ngenebio/pipeline/step_modules/DnaAlignment/MarkDuplicates.wdl
#================================================================================

workflow markDuplicatesRun {
    String inputBam
    String logFile
    String sampleName
    String tmp
    Int cpu
    Int ram

    call markDuplicates{
        input:
        inputBam = inputBam,
        logFile = logFile,
        tmpDir = tmp,
        cpu = cpu,
        ram = ram,
    }
}


task markDuplicates{
    ## Input Paramter
    String stepName
    String inputBam
    String logFile
    String tmpDir
    Int cpu
    Int ram

    # Set tools
    String java = '/opt/ngenebio/app/java/jre1.8.0_65/bin/java'
    String picard = '/opt/ngenebio/app/picard/picard-tools-2.5.0/picard.jar'

    # Set Parameter
    String outputDpBam = sub(inputBam, '\\.raw.bam$', '.dp.bam')
    String outputDpMetrix = sub(inputBam, '\\.raw.bam$', '.dp.metrix')
    
    command <<<
        export inputStep=${stepName}
        export logFile=${logFile}
        source /opt/ngenebio/pipeline/script/validate.sh

        ret1=$(md5sum_check ${outputDpBam})
        ret2=$(md5sum_check ${outputDpMetrix})
        if [[ $ret1 = 'True' && $ret2 = 'True' ]]; then
            log_progress 'Mark duplicates already finished!!!'
        else
            log_progress 'Mark duplicates start'
            wrap ${java} -XX:ParallelGCThreads=${cpu} -Xms${ram}g -Xmx${ram}g                              \
                -Dsamjdk.compression_level=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10                              \
                -jar ${picard} MarkDuplicates INPUT=${inputBam} OUTPUT=${outputDpBam} TMP_DIR=${tmpDir}   \
                METRICS_FILE=${outputDpMetrix} MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=200 CREATE_INDEX=true              \
                CREATE_MD5_FILE=false COMPRESSION_LEVEL=1 VALIDATION_STRINGENCY=SILENT

            wrap md5sum ${outputDpBam} > '${outputDpBam}.md5'
            wrap md5sum ${outputDpMetrix} > '${outputDpMetrix}.md5'
            log_progress 'Mark duplicates finished'
        fi
    >>>

    output {
        String dpBamFile = outputDpBam
        String dpMetrixFile = outputDpMetrix
    }

}