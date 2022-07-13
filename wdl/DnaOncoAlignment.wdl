#================================================================================
#
# FILE: DnaOncoAlignment.wdl
#
# DESCRIPTION: NGeneBio DNA ONCOaccuPanel DnaOncoAlignment workflow
# Note:
# REQUIREMENTS:
# BUGS:
# NOTES:
# AUTHOR: Juyeon hong hong.juyeon@ngenebio.com
# COMPANY: ngenebio
# VERSION: 0.1
# CREATED: 2022.04.24
# REVISION:
# RUN COMMAND: /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/cromwell-50.jar run /opt/ngenebio/pipeline/wdl/DnaOncoAlignment.wdl -i /opt/ngenebio/test/testOncoAlignment.json
# VALIDATION : /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/womtool-50.jar validate /opt/ngenebio/pipeline/wdl/DnaOncoAlignment.wdl -i /opt/ngenebio/test/testOncoAlignment.json
#================================================================================

import '/opt/ngenebio/pipeline/step_modules/DnaAlignment/BwaV059.wdl' as BwaV059
import '/opt/ngenebio/pipeline/step_modules/DnaAlignment/MarkDuplicates.wdl' as MarkDup
import '/opt/ngenebio/pipeline/step_modules/DnaAlignment/GATKV1650.wdl' as GATKV1650
import '/opt/ngenebio/pipeline/step_modules/DnaAlignment/CopyFinalBam.wdl' as CopyFinalBam
import '/opt/ngenebio/pipeline/wdl/Setting.wdl' as Setting

workflow DnaOncoAlignment {

    # Input Parameter
    String stepName='Alignment'
    String inputFasta
    String inputFastqR1 ## copied FastqR1 (output of DnaFastqValidation)
    String inputFastqR2 ## copied FastqR2 (output of DnaFastqValidation)
    String sampleName
    Int sampleId
    Int runId

    String fastqDir #/opt/ngenebio/output/runId/fastq/
    String logDir #/opt/ngenebio/output/runId/logs/
    String resultDir
    String alignDir
    String tmpDir

    String inputIntervalBaits  ## Workflow에서 input으로 받음 version에 따라서 달라짐
    String inputMills = '/opt/ngenebio/dependencies/1KGP/Mills_and_1000G_gold_standard.indels.b37.ngb.vcf.gz'
    String inputG1000Phase1 = '/opt/ngenebio/dependencies/1KGP/1000G_phase1.indels.b37.ngb.vcf.gz'
    String dbsnp137 = '/opt/ngenebio/dependencies/dbSNP/dbsnp_137.b37.vcf.gz'
    Int cpu = 2
    Int ram = 4

    String logFile = logDir + sampleId + '.' + runId + '.log.align.stderr.txt'
    String finalOutput  = resultDir + sampleName + '.final.bam'

    call Setting.checkMd5Sum {
        input:
            logFile = logFile,
            inputFile = finalOutput,
            stepName = stepName
    }

    if (checkMd5Sum.workMd5Sum == 'true') {
        call Setting.pseudoOutput {
            input:
                inputFile = finalOutput
        }
    }

    if (checkMd5Sum.workMd5Sum != 'true') {

        # FastqR1
        call BwaV059.bwaAln as runBwaAlnR1 {
            input:
                stepName = stepName + '-bwaAln-R1',
                readStrand = 'R1',
                logFile = logFile,
                inputFasta = inputFasta,
                outputDir = alignDir,
                inputFastq = inputFastqR1,
                cpu = cpu
        }

        # FastqR2
        call BwaV059.bwaAln as runBwaAlnR2 {
            input:
                stepName = stepName + '-bwaAln-R2',
                readStrand = 'R2',
                logFile = logFile,
                inputFasta = inputFasta,
                outputDir = alignDir,
                inputFastq = inputFastqR2,
                cpu = cpu
        }

        # Raw Bam File
        call BwaV059.bwaSampe as runSampe {
            input:
                stepName = stepName + '-bwaSampe',
                logFile = logFile,
                inputFastqR1 = inputFastqR1,
                inputFastqR2 = inputFastqR2,
                inputSaiR1 = runBwaAlnR1.saiFile,
                inputSaiR2 = runBwaAlnR2.saiFile,
                inputFasta = inputFasta,
                outputDir = alignDir,
                sampleName = sampleName,
                tmpDir = tmpDir,
                cpu = cpu,
                ram = ram
        }

        # MarkDuplicates
        call MarkDup.markDuplicates as markDuplicates {
            input:
                stepName = stepName + '-markDuplicates',
                inputBam = runSampe.bamFile,
                logFile = logFile,
                tmpDir = tmpDir,
                cpu = cpu,
                ram = ram
        }

        # IndelRealign
        call GATKV1650.indelRealign as realignIndel {
            input:
                stepName = stepName + '-indelRealign',
                inputBam = markDuplicates.dpBamFile,
                inputFasta = inputFasta,
                inputIntervalBaits = inputIntervalBaits,
                inputMills = inputMills,
                inputG1000Phase1 = inputG1000Phase1,
                logFile = logFile,
                tmpDir = tmpDir,
                cpu = cpu,
                ram = ram
        }

        # CountCovariates
        call GATKV1650.countCovariates as countCovariates {
            input:
                stepName = stepName + '-countCovariate',
                inputBam = realignIndel.realignBamFile,
                inputFasta = inputFasta,
                dbsnp = dbsnp137,
                logFile = logFile,
                tmpDir = tmpDir,
                cpu = cpu,
                ram = ram
        }

        # BaseRecalibration
        call GATKV1650.baseRecalibration as recalBase {
            input:
                stepName = stepName + '-baseRecalibration',
                inputBam = realignIndel.realignBamFile,
                inputCsv = countCovariates.covariateCsvFile,
                inputFasta = inputFasta,
                logFile = logFile,
                tmpDir = tmpDir,
                cpu = cpu,
                ram = ram
        }

        call CopyFinalBam.copyFinalBam as copyFinalBam {
            input :
                stepName = stepName + '-copyFinalFile',
                inputBam = recalBase.recalBamFile,
                outputBam = finalOutput,
                logFile = logFile
        }
    }

    # Outputs that will be retained when execution is complete
    output {
            String finalBamFile = select_first([copyFinalBam.copyOutput, pseudoOutput.pseudoFile])
        }
}
