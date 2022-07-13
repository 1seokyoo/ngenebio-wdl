#================================================================================
#
# FILE: BrcaComplementSoftclip.wdl
#
# DESCRIPTION: Complement soft-clip using inhouse script (clip_finder)
# Note:
# REQUIREMENTS:
# BUGS:
# NOTES:
# AUTHOR: Wonseok Yoo wonseok.yoo@ngenebio.com
# COMPANY: ngenebio
# VERSION: 0.1
# CREATED: 2022.05.29
# REVISION:
# RUN COMMAND:
# VALIDATION :
#================================================================================

task findSoftclip {

    # Input Parameter
    String stepName
    String logFile
    String sampleName
    String inputFile
    String outputDir
    String outputName
    String primerBed

    String tmpSamFile = outputDir + outputName + '.noHdr.sam'
    String outputBedFile = outputDir + outputName + '.bed'
    String outputTxtFile = outputDir + outputName + '.txt'

    # Used tools
    String inhouseScript = '/opt/ngenebio/pipeline/step_modules/DnaAlignment/brca_findSoftclip.py'
    String samtools = '/opt/ngenebio/app/samtools/samtools-1.3.1/samtools'

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

        ret1=$(md5sum_check ${outputBedFile})
        ret2=$(md5sum_check ${outputTxtFile})
        if [[ $ret1 = 'True' && $ret2 = 'True' ]]; then
            log_progress 'Already finished to find soft-clip using inhouse script (clip_finder)'
        else
            log_progress 'Start to find soft-clip using inhouse script (clip_finder)'

            wrap ${samtools} view -F 258 ${inputFile} > ${tmpSamFile}
            wrap python ${inhouseScript} amplicon=${primerBed} \
                                         input_sam=${tmpSamFile} \
                                         sample_name=${sampleName} \
                                         output_bed=${outputBedFile} \
                                         output_txt=${outputTxtFile} \
                                         1>&2

            wrap md5sum ${outputBedFile} > '${outputBedFile}.md5'
            wrap md5sum ${outputTxtFile} > '${outputTxtFile}.md5'
            log_progress 'Finish to find soft-clip using inhouse script (clip_finder)'
        fi

        outputTxt='false'
        if [ -s ${outputBedFile} ]; then
            log_progress "BWA realign with warnning bed file"
            outputTxt='true'
        fi

        echo $outputTxt
    >>>

    output {
        String bedExist = read_string(stdout())
        String bedFile = outputBedFile
        String txtFile = outputTxtFile
    }
}

task extractSoftclipFastq {

    # Input Parameter
    String stepName
    String logFile
    String sampleName
    String inputFile
    String inputBed
    String outputDir
    String outputName

    String outputBamExtractFile = outputDir + outputName + '.extract.bam'
    String outputBamRemainFile = outputDir + outputName + '.remain.bam'
    String outputFastqExtractFile = outputDir + outputName + '.extract.fastq'

    # Used tools
    String inhouseScript = '/opt/ngenebio/pipeline/step_modules/DnaAlignment/brca_splitBam.py'
    String bedtools = '/opt/ngenebio/app/bedtools/bedtools-v2.26.0/bin/bedtools'

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

        ret=$(md5sum_check ${outputFastqExtractFile})
        if [[ $ret = 'True' ]]; then
            log_progress 'Already finished to extract fastq in soft-clip region'
        else
            log_progress 'Start to extract fastq in soft-clip region'

            wrap python ${inhouseScript} ${inputBed} \
                                         ${inputFile} \
                                         ${outputBamExtractFile} \
                                         ${outputBamRemainFile}

            wrap ${bedtools} bamtofastq -i ${outputBamExtractFile} \
                                        -fq ${outputFastqExtractFile}

            wrap md5sum ${outputFastqExtractFile} > '${outputFastqExtractFile}.md5'
            log_progress 'Finish to extract fastq in soft-clip region'
        fi
    >>>

    output {
        String fastqFile = outputFastqExtractFile
        String retainFile = outputBamRemainFile
    }
}

task mergeBam {

    # Input Parameter
    String stepName
    String logFile
    String inputFile1
    String inputFile2
    String outputDir
    String outputName
    String tmpDir

    String outputRawFile = outputDir + outputName + '.merged.raw.bam'
    String outputSortFile = outputDir + outputName + '.merged.sort.bam'

    # Used tools
    String samtools = '/opt/ngenebio/app/samtools/samtools-1.3.1/samtools'

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

        ret=$(md5sum_check ${outputSortFile})
        if [[ $ret = 'True' ]]; then
            log_progress 'Already finished to merge soft-clip bam files'
        else
            log_progress 'Start to merge soft-clip bam files'

            wrap ${samtools} merge -f ${outputRawFile} \
                                      ${inputFile1} \
                                      ${inputFile2}

            wrap ${samtools} sort -O bam -T ${tmpDir} \
                                  -o ${outputSortFile} \
                                  ${outputRawFile}


            wrap ${samtools} index ${outputSortFile}

            wrap md5sum ${outputSortFile} > '${outputSortFile}.md5'
            log_progress 'Finish to merge soft-clip bam files'
        fi
    >>>

    output {
        String bamFile = outputSortFile
    }
}

task runGmap {

    # Input Parameter
    String stepName
    String logFile
    String sampleName
    String inputFile
    String outputDir
    String tmpDir
    String cpu

    String outputRawSamFile = outputDir + sampleName + '.gmap.raw.sam'
    String outputDelSamFile = outputDir + sampleName + '.gmap.sam'

    # Used tools
    String gmap = '/opt/ngenebio/app/gmap-2016-11-07/bin/gmap'
    String gmap_db = '/opt/ngenebio/dependencies/gmap'

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

        ret=$(md5sum_check ${outputDelSamFile})
        if [[ $ret = 'True' ]]; then
            log_progress 'Already finished to align using GMAP'
        else
            log_progress 'Start to align using GMAP'

            wrap ${gmap} -t ${cpu} \
                         -D ${gmap_db} -d brca \
                         -f samse -n 0 \
                         ${inputFile} 2> ${outputRawSamFile}

            awk '{gsub("N","D",$6)}1' ${outputRawSamFile} | awk -v OFS="\t" '$1=$1' > ${outputDelSamFile}

            wrap md5sum ${outputDelSamFile} > '${outputDelSamFile}.md5'
            log_progress 'Finish to align using GMAP'
        fi
    >>>

    output {
        String samFile = outputDelSamFile
    }
}