#================================================================================
#
# FILE: Setting.wdl
#
# DESCRIPTION: NGeneBio setting task
# Note:
# REQUIREMENTS:
# BUGS:
# NOTES:
# AUTHOR: Wonseok Yoo wonseok.yoo@ngenebio.com
# COMPANY: ngenebio
# VERSION: 1.0
# CREATED: 2022.04.24
# REVISION:
# VALIDATION: /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/womtool-50.jar validate /opt/ngenebio/pipeline/wdl/Setting.wdl
#================================================================================

task setDir {

    # Input Parameter
    String runId

    # Set Parameter
    String runDir = '/opt/ngenebio/output/' + runId + '/'
    String fastqDir = runDir + 'fastq/'
    String dataDir = runDir + 'data/'
    String logsDir = runDir + 'logs/'
    String tmpDir = runDir + 'tmp/'
    String assayDir = runDir + 'assay_reference/'
    String resultDir = runDir + 'result/'
    String fastqcDir = dataDir + 'fastqc/'
    String alignDir = dataDir + 'alignment/'
    String statDir = dataDir + 'stat/'
    String snvDir = dataDir + 'variant/snv/'
    String cnvtDir = dataDir + 'variant/cnv/'
    String svDir = dataDir + 'variant/sv/'
    
    command <<<
        
        # Ensure script halts when non-zero exit status detected
        set -e
        set -o pipefail # trace ERR through pipes
        set -o errtrace  # trace ERR through 'time command' and other functions
        set -o nounset   # set -u : exit the script if you try to use an uninitialised variable
        set -o errexit   # set -e : exit the script if any statement returns a non-true return value
        ulimit -n 2048 # Open file count

        mkdir -p ${fastqDir}
        mkdir -p ${dataDir}
        mkdir -p ${logsDir}
        mkdir -p ${tmpDir}
        mkdir -p ${snvDir}
        mkdir -p ${cnvtDir}
        mkdir -p ${svDir}
        mkdir -p ${resultDir}
        mkdir -p ${fastqcDir}
        mkdir -p ${alignDir}
        mkdir -p ${statDir}
        mkdir -p ${assayDir}
        rm -r -f ${logsDir}/*
        
    >>>

    output {
        String run = runDir
        String fastq = fastqDir
        String fastqc = fastqcDir
        String data = dataDir
        String logs = logsDir
        String tmp = tmpDir
        String snv = snvDir
        String cnv = cnvtDir
        String sv = svDir
        String result = resultDir
        String stat = statDir
        String align = alignDir
        String assay = assayDir
    }
}


task checkMd5Sum {

    String logFile
    String stepName
    String inputFile

    command <<<
        export inputStep=${stepName}
        export logFile=${logFile}
        source /opt/ngenebio/pipeline/script/validate.sh
        
        ret=$(md5sum_check ${inputFile})
        workMd5Sum='false'

        if [[ $ret = 'True' ]]; then

            log_progress "${stepName} already finished !!!"

            workMd5Sum='true'
            
        fi

        echo $workMd5Sum
        
    >>>

    output {
        String workMd5Sum = read_string(stdout())
    }

}


task pseudoOutput {

    String inputFile
    String outputFile = inputFile

    command <<<
    >>>

    output {
        String pseudoFile = outputFile
    }

}


task splitPipelineCode {
    
    String pipelineCode 

    command <<<
        python <<CODE
        
        pipelineInfoList = '${pipelineCode}'.split("-")
        idx = range(0,len(pipelineInfoList))
        for i in idx :
            print(pipelineInfoList[i])

        CODE
    >>>
    output {
        Array[String] pipelineInfo = read_lines(stdout())
    }
}


task setBrcaTarget {

    String outputDir
    String assayTargetDir

    String panelAmpliconBed = assayTargetDir + 'amplicon.bed'
    String panelPrimerBed = assayTargetDir + 'primer.bed'
    String panelPrimerBedpe = assayTargetDir + 'primer.bedpe'
    String panelRoiBed = assayTargetDir + 'roi.bed'

    String sampleAmpliconBed = outputDir + 'amplicon.bed'
    String samplePrimerBed = outputDir + 'primer.bed'
    String samplePrimerBedpe = outputDir + 'primer.bedpe'
    String sampleRoiBed = outputDir + 'roi.bed'

    command <<<

        # Ensure script halts when non-zero exit status detected
        set -e
        set -o pipefail # trace ERR through pipes
        set -o errtrace  # trace ERR through 'time command' and other functions
        set -o nounset   # set -u : exit the script if you try to use an uninitialised variable
        set -o errexit   # set -e : exit the script if any statement returns a non-true return value
        ulimit -n 2048 # Open file count

        cp ${assayTargetDir}/* ${outputDir}
        
    >>>

    output {
        String ampliconBed = sampleAmpliconBed
        String primerBed = samplePrimerBed
        String primerBedpe = samplePrimerBedpe
        String roiBed = sampleRoiBed
    }
}


task setOncoTarget {

    String outputDir
    String assayDir
    String pipelineVersion

    String assayTargetDir = assayDir + pipelineVersion

    String panelHrdGeneList = assayDir + 'COMMON/hrdGeneList.txt'
    String panelMmrGeneList = assayDir + 'COMMON/mmrGeneList.txt'
    String panelCnvReference = '/opt/ngenebio/dependencies/CNVkit/ONCO/' + pipelineVersion + '/AMC-NORMAL-00001-00004.cnn'
    String panelSvReference =  '/opt/ngenebio/dependencies/BreaKmer/ONCO/' + pipelineVersion + '/' 

    String sampleBaitsFile = outputDir + '/baitsInterval.list'
    String sampleSnvAllBed = outputDir + '/snvAllGene.bed'
    String sampleSnvHrdBed = outputDir + '/snvHrdGene.bed'
    String sampleSnvWhiBed = outputDir + '/snvWhiGene.bed'
    String sampleSnvWhiTxt = outputDir + '/snvWhiGene.txt'
    String sampleSnvTmpVcf = outputDir + '/snvTmp.vcf'
    String sampleSvTagetBed = outputDir + '/svTarget.bed'
    String sampleSvTagetTxt = outputDir + '/svTarget.txt'
    String sampleHrdGeneList = outputDir + '/hrdGeneList.txt'
    String sampleMmrGeneList = outputDir + '/mmrGeneList.txt'

    command <<<

        # Ensure script halts when non-zero exit status detected
        set -e
        set -o pipefail # trace ERR through pipes
        set -o errtrace  # trace ERR through 'time command' and other functions
        set -o nounset   # set -u : exit the script if you try to use an uninitialised variable
        set -o errexit   # set -e : exit the script if any statement returns a non-true return value
        ulimit -n 2048 # Open file count

        cp ${assayTargetDir}/* ${outputDir}
        cp ${panelHrdGeneList} ${sampleHrdGeneList}
        cp ${panelMmrGeneList} ${sampleMmrGeneList}
        
    >>>

    output {
        String baitsFile = sampleBaitsFile
        String snvAllBed = sampleSnvAllBed
        String snvHrdBed = sampleSnvHrdBed
        String snvWhiBed = sampleSnvWhiBed
        String snvWhiTxt = sampleSnvWhiTxt
        String snvTmpVcf = sampleSnvTmpVcf
        String svTagetBed = sampleSvTagetBed
        String svTagetTxt = sampleSvTagetTxt
        String hrdGeneList = sampleHrdGeneList
        String mmrGeneList = sampleMmrGeneList
        String cnvReference = panelCnvReference
        String svReference = panelSvReference
    }
}


task copyFinalFile {

    String stepName
    String logFile
    String inputFile
    String outputName

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

        ret=$(md5sum_check ${outputName})

        if [[ $ret = 'True' ]]; then
            log_progress 'Final File already copied to results Directory !!!'

        else
            log_progress 'Copy Final File to results Directory'

            cp ${inputFile} ${outputName}
            md5sum ${outputName} > '${outputName}.md5'

            log_progress 'Copy Final File to results Directory Done'
        fi
    >>>

    output {
        String outputFile = outputName
    }
}
