#================================================================================
#
# FILE: BRCAPreprocessing.wdl
#
# DESCRIPTION: NGeneBio DNA BRCAaccuPanel Fastq Preprocessing workflow
# Note:
# REQUIREMENTS:
# BUGS:
# NOTES:
# AUTHOR: Wonseok Yoo wonseok.yoo@ngenebio.com
# COMPANY: ngenebio
# VERSION: 0.1
# CREATED: 2022.05.24
# REVISION:
# RUN COMMAND: 
# VALIDATION : 
#================================================================================

task normalizeReadDepth {

    # Input parameter
    String logFile
    String stepName
    String outputDir
    String sampleName
    String normalizeDir = outputDir + sampleName + '_normalize'

    # Fastq file name
    String inputFastqR1
    String inputFastqR2
    String tmpFastqR1 = sub(inputFastqR1, "\\.fastq.gz", "\\.tmp.fastq")
    String tmpFastqR2 = sub(inputFastqR2, "\\.fastq.gz", "\\.tmp.fastq")
    String outputFastqR1 = outputDir + sampleName + '_normalize.1.fastq'
    String outputFastqR2 = outputDir + sampleName + '_normalize.2.fastq'
    
    # Normalize setting
    Int cycle
    Int expectedDepth
    Int targetLength

    # Used tools
    String fastqToolsSample = '/opt/ngenebio/app/fastq-tools-0.8/fastq-sample'

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

        ret1=$(md5sum_check ${outputFastqR1})
        ret2=$(md5sum_check ${outputFastqR2})
        if [[ $ret1 = 'True' && $ret2 = 'True' ]]; then
            log_progress 'Already finished to normalize fastq'
        else
            log_progress 'Start to normalize fastq'

            expectedFqLine=`expr ${targetLength} / ${cycle} \* ${expectedDepth}`
            expectedPairFqLine=`expr $expectedFqLine / 2`

            wrap gzip -dc ${inputFastqR1} > ${tmpFastqR1}
            wrap gzip -dc ${inputFastqR2} > ${tmpFastqR2}

            wrap ${fastqToolsSample} -n $expectedPairFqLine ${tmpFastqR1} ${tmpFastqR2} -o ${normalizeDir}
            
            wrap rm -r -f ${tmpFastqR1}
            wrap rm -r -f ${tmpFastqR2}
            wrap md5sum ${outputFastqR1} > '${outputFastqR1}.md5'
            wrap md5sum ${outputFastqR2} > '${outputFastqR2}.md5'
            
            log_progress 'Finish to normalize fastq'
        fi

    >>>

    output {
        String fastqR1File = outputFastqR1
        String fastqR2File = outputFastqR2
    }
}


task removePrimer {

    # Input parameter
    String logFile
    String stepName
    String outputDir 
    String sampleName
    String primerBed

    String inputFastqR1
    String inputFastqR2

    String assayDir
    String ampBed = assayDir + sampleName + '_cut_amp.fasta'
    String cutPrimer = assayDir + sampleName + '_cutted_amp.fasta'

    # Output Parameter
    String fastqSampleName = outputDir + sampleName
    String ahoOut1 = fastqSampleName + '_1.aho'
    String ahoOut2 = fastqSampleName + '_2.aho'
    String unknownOut = fastqSampleName + '_cut_5bp_primer'
    String unknownOut2 = fastqSampleName + '_un_2'
    String unknownOut3 = fastqSampleName + '_un_3'
    String priStat = fastqSampleName + '_primer.stat'
    String ampliconStatOut = fastqSampleName + '_fastq_ampstat.txt'
    String perfectReads = fastqSampleName + '_perfect.list'
    String outputFastqR1 = fastqSampleName + '_normalize.1_trim.fastq'
    String outputFastqR2 = fastqSampleName + '_normalize.2_trim.fastq'
    String cutadaptFastqR1 = fastqSampleName + '_normalize.1_cutadapt.fastq' #fastq 에서 fasta 로 바꿀 수 있으면 변경
    String cutadaptFastqR2 = fastqSampleName + '_normalize.2_cutadapt.fastq'
    String notoutputFastqR1 = fastqSampleName + '_normalize.1_untrimmed.fastq'
    String notoutputFastqR2 = fastqSampleName + '_normalize.2_untrimmed.fastq'

    # Used tools
    String java = "/opt/ngenebio/app/java/jre1.8.0_144/bin/java"
    String bedToBed = "/opt/ngenebio/app/primer_remover/bed_to_bed.py" #bedtobed -> primerremover mv
    String primerNtCut = "/opt/ngenebio/app/primer_remover/primer_nt_cut.py"
    String exactRemover = "/opt/ngenebio/app/primer_remover/Exact_remover.jar"
    String primerInput2 = "/opt/ngenebio/app/primer_remover/primer_input2.py"
    String newFastq = "/opt/ngenebio/app/primer_remover/primer_remove_fastq_NGB_SW.py"
    String ampliconCheckV3 = "/opt/ngenebio/app/primer_remover/amplicon_check_v3.py"
    String filter = "/opt/ngenebio/app/primer_remover/filter.py"

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

        ret1=$(md5sum_check ${outputFastqR1})
        ret2=$(md5sum_check ${outputFastqR2})
        if [[ $ret1 = 'True' && $ret2 = 'True' ]]; then
            log_progress 'Already finished to find perfect fastq'
        else
            log_progress 'Start to find perfect fastq'

            wrap python ${bedToBed} bed_file=${primerBed} bed_type=primer_cut > ${ampBed}
            wrap python ${primerNtCut} fastq=${ampBed} cut_num=5 output=${cutPrimer}
            wrap ${java} -jar ${exactRemover} ${inputFastqR1} ${cutPrimer} > ${ahoOut1}
            wrap ${java} -jar ${exactRemover} ${inputFastqR2} ${cutPrimer} > ${ahoOut2}
            wrap python ${primerInput2} cut_primer=${cutPrimer} output=${unknownOut}
            (wrap python ${newFastq} ${ahoOut1} ${inputFastqR1} ${unknownOut} ${primerBed} > ${unknownOut2}) &
            (wrap python ${newFastq} ${ahoOut2} ${inputFastqR2} ${unknownOut} ${primerBed} > ${unknownOut3}) &
            wait
            wrap python ${ampliconCheckV3} ${unknownOut2} ${unknownOut3} ${primerBed} ${ampliconStatOut} > ${priStat}
            wrap cat ${priStat} | grep '^@' > ${perfectReads}
            (wrap python ${filter} ${perfectReads} ${outputFastqR1} ${cutadaptFastqR1} ${notoutputFastqR1}) &
            (wrap python ${filter} ${perfectReads} ${outputFastqR2} ${cutadaptFastqR2} ${notoutputFastqR2}) &
            wait
            wrap md5sum ${outputFastqR1} > '${outputFastqR1}.md5'
            wrap md5sum ${outputFastqR2} > '${outputFastqR2}.md5'
            log_progress 'Finish to find perfect fastq'
        fi
    >>>

    output {
        String fastqR1File = outputFastqR1
        String fastqR2File = outputFastqR2
        String cutadaptR1File = cutadaptFastqR1
        String cutadaptR2File = cutadaptFastqR2
        String untrimR1File = notoutputFastqR1
        String untrimR2File = notoutputFastqR2
    }
}