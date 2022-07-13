#================================================================================
#
# FILE: VardictV16Run.wdl
#
# DESCRIPTION: NGeneBio DNA ONCOaccuPanel Vardict variant calling workflow
# Note:
# REQUIREMENTS:
# BUGS:
# NOTES:
# AUTHOR: Juyeon hong hong.juyeon@ngenebio.com
# COMPANY: ngenebio
# VERSION: 0.1
# CREATED: 2022.04.24
# REVISION:
# RUN COMMAND:
# VALIDATION : /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/womtool-50.jar validate /opt/ngenebio/pipeline/step_modules/SnvCalling/VardictV16Run.wdl
#================================================================================

task vardictCall{

    # Input Parameter
    String stepName
    String inputTarget
    String inputBam
    String sampleName
    String inputFasta
    String mode 
    String outputDir
    String optionFreq
    Boolean optionGenomicCall  ## true or false
    String logFile
    Int cpu

    String outputVardict = outputDir + sampleName + '.' + mode + '.vardict'
    String genomicOption = if optionGenomicCall then "-p" else ""

    # Used tools
    String vardict = '/opt/ngenebio/app/VarDict-1.6/vardict'

    command <<<
        export inputStep=${stepName}
        export logFile=${logFile}
        source /opt/ngenebio/pipeline/script/validate.sh

        ret=$(md5sum_check ${outputVardict})
        if [[ $ret = 'True' ]]; then
            log_progress 'Variant Call (VarDict) already finished!!!'
        else
            log_progress 'Variant Call (VarDict) start'

            wrap source /snv-pipeline/bin/activate

            wrap ${vardict} -G ${inputFasta} -t ${genomicOption} -N ${sampleName} -b ${inputBam} -c 1 -S 2 -E 3 -g 5 -f ${optionFreq} -th ${cpu} ${inputTarget} > ${outputVardict}

            wrap md5sum ${outputVardict} > '${outputVardict}.md5'
            log_progress 'Variant Call (VarDict) finished'
        fi

    >>>

    output {
        String rawVardictFile = outputVardict
    }
}

task teststrandbias{

    # Input Parameter
    String stepName
    String inputVardict
    String outputDir
    String sampleName
    String mode
    String logFile
    Int cpu

    String outputStrandbias = outputDir + sampleName + '.' + mode + '.strandbias'

    # Used tools
    String script = '/opt/ngenebio/app/VarDict-1.6/teststrandbias.R'

    command <<<
        export inputStep=${stepName}
        export logFile=${logFile}
        source /opt/ngenebio/pipeline/script/validate.sh

        ret=$(md5sum_check ${outputStrandbias})
        if [[ $ret = 'True' ]]; then
            log_progress 'Variant Call (Strand Bias) already finished!!!'
        else
            log_progress 'Variant Call (Strand Bias) start'

            wrap source /snv-pipeline/bin/activate
            wrap cat ${inputVardict} | Rscript ${script} > ${outputStrandbias}

            wrap md5sum ${outputStrandbias} > '${outputStrandbias}.md5'
            log_progress 'Variant Call (Strand Bias) finished'
        fi
    >>>

    output {
        String rawStrandbiasFile = outputStrandbias
    }
}


task var2Vcf{

    # Input Parameter
    String stepName
    String inputStrandbias
    Int optionAo
    String optionFreq 
    String outputDir
    String sampleName
    String mode
    String logFile
    Int cpu

    String outputVcf = outputDir + sampleName + '.' + mode + '.raw.vcf'

    # Used tools
    String var2Vcf = '/opt/ngenebio/app/VarDict-1.6/var2vcf_valid.pl'

    command <<<
        export inputStep=${stepName}
        export logFile=${logFile}
        source /opt/ngenebio/pipeline/script/validate.sh

        ret=$(md5sum_check ${outputVcf})
        if [[ $ret = 'True' ]]; then
            log_progress 'Variant Call (Vardic to VCF) already finished!!!'
        else
            log_progress 'Variant Call (Vardict to VCF) start'

            wrap source /snv-pipeline/bin/activate
            (wrap ${var2Vcf} -N ${sampleName} -Q 20 -d 30 -v ${optionAo} -f ${optionFreq} ${inputStrandbias} > ${outputVcf}) &
            wait

            wrap md5sum ${outputVcf} > '${outputVcf}.md5'
            log_progress 'Variant Call (Vardic to VCF) finished'
        fi

    >>>

    output {
        String rawVcfFile = outputVcf
    }
}


