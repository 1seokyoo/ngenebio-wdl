#================================================================================
#
# FILE: CopyFinalBam.wdl
#
# DESCRIPTION: NGeneBio DNA ONCOaccuPanel Final Bam Copy workflow
# Note:
# REQUIREMENTS:
# BUGS:
# NOTES:
# AUTHOR: Juyeon hong hong.juyeon@ngenebio.com
# COMPANY: ngenebio
# VERSION: 0.1
# CREATED: 2022.05.10
# REVISION:
# RUN COMMAND: /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/cromwell-50.jar run /opt/ngenebio/pipeline/step_modules/VariantCalling/MakeWhitelistTarget.wdl -i /opt/ngenebio/test/testONCOalignment.json
# VALIDATION : /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/womtool-50.jar validate /opt/ngenebio/pipeline/step_modules/VariantCalling/MakeWhitelistTarget.wdl -i /opt/ngenebio/test/testONCOalignment.json
#================================================================================

task copyFinalBam {

    String inputBam
    String stepName
    String outputBam 
    String inputBai = inputBam + '.bai'
    String outputBai = outputBam + '.bai'
    String logFile

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

        ret=$(md5sum_check ${outputBam})

        if [[ $ret = 'True' ]]; then
            log_progress 'Final Bam File already copied to results Directory !!!'

        else
            log_progress 'Copy Final Bam File to results Directory'

            cp ${inputBam} ${outputBam}
            cp ${inputBai} ${outputBai}

            wrap md5sum ${outputBam} > '${outputBam}.md5'

            log_progress 'Copy Final Bam File to results Directory Done'
        fi
    >>>

    output {

        String copyOutput = outputBam

    }

}