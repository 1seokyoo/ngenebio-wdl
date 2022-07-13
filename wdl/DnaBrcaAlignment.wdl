#================================================================================
#
# FILE: DnaBrcaAlignment.wdl
#
# DESCRIPTION: NGeneBio DNA BRCAaccuTest DnaBrcaAlignment workflow
# Note:
# REQUIREMENTS:
# BUGS:
# NOTES:
# AUTHOR: Wonseok Yoo wonseok.yoo@ngenebio.com
# COMPANY: ngenebio
# VERSION: 0.1
# CREATED: 2022.04.24
# REVISION:
# RUN COMMAND: /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/cromwell-50.jar run /opt/ngenebio/pipeline/wdl/DnaBrcaAlignment.wdl -i /opt/ngenebio/test/testBrcaAlignment.json
# VALIDATION : /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/womtool-50.jar validate /opt/ngenebio/pipeline/wdl/DnaBrcaAlignment.wdl -i /opt/ngenebio/test/testBrcaAlignment.json
#================================================================================

import '/opt/ngenebio/pipeline/step_modules/DnaAlignment/BwaV0717.wdl' as BwaV0717
import '/opt/ngenebio/pipeline/step_modules/DnaAlignment/BrcaComplementSoftclip.wdl' as ComplementSoftclip
import '/opt/ngenebio/pipeline/step_modules/DnaAlignment/GATKV4261.wdl' as GATKV4261
import '/opt/ngenebio/pipeline/step_modules/DnaAlignment/BamclipperV111.wdl' as BamclipperV111
import '/opt/ngenebio/pipeline/step_modules/DnaAlignment/CopyFinalBam.wdl' as CopyFinalBam
import '/opt/ngenebio/pipeline/wdl/Setting.wdl' as Setting

workflow DnaBrcaAlignment {

    # Input Parameter
    String stepName = 'Alignment'
    String inputFastq
    String sampleName
    Int sampleId
    Int runId
    String resultDir
    String alignDir #/opt/ngenebio/output/runId/data/alignment/
    String logDir #/opt/ngenebio/output/runId/logs/
    String tmpDir
    String primerBed
    String primerBedpe
    Int cpu = 2
    Int memory = 4

    String referenceFasta
    String inputMills = '/opt/ngenebio/dependencies/1KGP/Mills_and_1000G_gold_standard.indels.b37.ngb.vcf.gz'
    String inputG1000Phase1 = '/opt/ngenebio/dependencies/1KGP/1000G_phase1.indels.b37.ngb.vcf.gz'
    String dbsnp155 = '/opt/ngenebio/dependencies/dbSNP/dbsnp_155_all.vcf.gz'

    String logFile = logDir + sampleId + '.' + runId + '.log.align.stderr.txt'
    String bqsrOutput  = alignDir + sampleName + '.bqsr.bam'
    String finalOutput  = resultDir + sampleName + '.final.bam'

    call Setting.checkMd5Sum {
        input:
            stepName = stepName + '-checkMd5Sum',
            inputFile = finalOutput,
            logFile = logFile
    }

    if (checkMd5Sum.workMd5Sum == 'true') {
        call Setting.pseudoOutput as pseudoFinal {
            input:
                inputFile = finalOutput
        }
        call Setting.pseudoOutput as pseudoBqsr {
            input:
                inputFile = bqsrOutput
        }
    }

    if (checkMd5Sum.workMd5Sum != 'true') {

        call BwaV0717.bwaMem as runBwaMemDefault {
            input:
                stepName = stepName + '-bwaMemDefault',
                logFile = logFile,
                inputFile = inputFastq,
                outputDir = alignDir,
                outputName = sampleName,
                tmpDir = tmpDir,
                referenceFasta = referenceFasta,
                cpu = cpu
        }

        call ComplementSoftclip.findSoftclip as findSoftclipDefault {
            input:
                stepName = stepName + '-findSoftclipDefault',
                logFile = logFile,
                sampleName = sampleName,
                inputFile = runBwaMemDefault.samFile,
                outputDir = alignDir,
                outputName = sampleName + '.default',
                primerBed = primerBed
        }

        if (findSoftclipDefault.bedExist != 'true') {
            call GATKV4261.addOrReplaceReadGroups as sortBamDefault {
                input:
                    stepName = stepName + '-sortBamDefault',
                    logFile = logFile,
                    sampleName = sampleName,
                    inputFile = runBwaMemDefault.samFile,
                    outputDir = alignDir,
                    outputName = sampleName,
                    tmpDir = tmpDir,
                    cpu = cpu,
                    memory = memory    
            }
        }

        if (findSoftclipDefault.bedExist == 'true') {
            call ComplementSoftclip.extractSoftclipFastq as extractFastqFirst {
                input:
                    stepName = stepName + '-extractFastqFirst',
                    logFile = logFile,
                    sampleName = sampleName,
                    inputFile = runBwaMemDefault.samFile,
                    inputBed = findSoftclipDefault.bedFile,
                    outputDir = alignDir,
                    outputName = sampleName + '.extractFirst'
            }

            call BwaV0717.bwaMem as runBwaMemExtractFastq {
                input:
                    stepName = stepName + '-bwaMemExtractFastq',
                    logFile = logFile,
                    inputFile = extractFastqFirst.fastqFile,
                    outputDir = alignDir,
                    outputName = sampleName + '.extractFirst',
                    tmpDir = tmpDir,
                    referenceFasta = referenceFasta,
                    cpu = cpu,
                    mismatch_override = 2,
                    gap_override = 2
            }

            call ComplementSoftclip.findSoftclip as findSoftclipExtractFirst {
                input:
                    stepName = stepName + '-findSoftclipExtractFastq',
                    logFile = logFile,
                    sampleName = sampleName,
                    inputFile = runBwaMemExtractFastq.samFile,
                    outputDir = alignDir,
                    outputName = sampleName + '.extractFirst',
                    primerBed = primerBed
            }

            if (findSoftclipExtractFirst.bedExist != 'true') {
                call GATKV4261.addOrReplaceReadGroups as sortBamExtractFirst {
                    input:
                        stepName = stepName + '-sortBamExtractFirst',
                        logFile = logFile,
                        sampleName = sampleName,
                        inputFile = runBwaMemDefault.samFile,
                        outputDir = alignDir,
                        outputName = sampleName + '.extractFirst',
                        tmpDir = tmpDir,
                        cpu = cpu,
                        memory = memory    
                }

                call ComplementSoftclip.mergeBam as mergeBamExtractFirst {
                    input:
                        stepName = stepName + '-mergeBamExtractFirst',
                        logFile = logFile,
                        inputFile1 = extractFastqFirst.retainFile,
                        inputFile2 = sortBamExtractFirst.bamFile,
                        outputDir = alignDir,
                        outputName = sampleName + '.extractFirst',
                        tmpDir = tmpDir
                }
            }

            if (findSoftclipExtractFirst.bedExist == 'true') {
                call ComplementSoftclip.extractSoftclipFastq as extractFastqSecond {
                    input:
                        stepName = stepName + '-extractFastqSecond',
                        logFile = logFile,
                        sampleName = sampleName,
                        inputFile = runBwaMemDefault.samFile,
                        inputBed = findSoftclipExtractFirst.bedFile,
                        outputDir = alignDir,
                        outputName = sampleName + '.extractSecond'
                }

                call ComplementSoftclip.runGmap {
                    input:
                        stepName = stepName + '-extractGmap',
                        logFile = logFile,
                        sampleName = sampleName,
                        inputFile = extractFastqSecond.fastqFile,
                        outputDir = alignDir,
                        tmpDir = tmpDir,
                        cpu = cpu
                }

                call GATKV4261.addOrReplaceReadGroups as sortBamGmap {
                    input:
                        stepName = stepName + '-sortBamGmap',
                        logFile = logFile,
                        sampleName = sampleName,
                        inputFile = runGmap.samFile,
                        outputDir = alignDir,
                        outputName = sampleName + '.extractGmap',
                        tmpDir = tmpDir,
                        cpu = cpu,
                        memory = memory  
                }

                call ComplementSoftclip.mergeBam as mergeBamExtractSecond {
                    input:
                        stepName = stepName + '-mergeBamExtractSecond',
                        logFile = logFile,
                        inputFile1 = extractFastqFirst.retainFile,
                        inputFile2 = sortBamGmap.bamFile,
                        outputDir = alignDir,
                        outputName = sampleName + '.extractSecond',
                        tmpDir = tmpDir
                }
            }
        }

        String rawBamFile = select_first([sortBamDefault.bamFile, mergeBamExtractFirst.bamFile, mergeBamExtractSecond.bamFile])
        # 1순위 sortBamDefault.bamFile
        # 2순위 mergeBamExtractFirst.bamFile
        # 3순위 mergeBamExtractSecond.bamFile

        call GATKV4261.baseRecalibration {
            input:
                stepName = stepName + '-baseRecalibration',
                logFile = logFile,
                sampleName = sampleName,
                inputFile = rawBamFile,
                outputDir = alignDir,
                tmpDir = tmpDir,
                referenceFasta = referenceFasta,
                dbsnp = dbsnp155,
                mills = inputMills,
                g1000Phase1 = inputG1000Phase1,
                cpu = cpu,
                memory = memory
        }


        call GATKV4261.applyBqsr {
            input:
                stepName = stepName + '-applyBQSR',
                logFile = logFile,
                sampleName = sampleName,
                inputBamFile = rawBamFile,
                inputRecal = baseRecalibration.recalFile,
                outputDir = alignDir,
                tmpDir = tmpDir,
                referenceFasta = referenceFasta,
                cpu = cpu,
                memory = memory
        }


        call BamclipperV111.bamclipper {
            input:
                stepName = stepName + '-bamclipper',
                logFile = logFile,
                sampleName = sampleName,
                inputFile = applyBqsr.bamFile,
                outputDir = alignDir,
                primerBedpe = primerBedpe,
                cpu = cpu
        }


        call CopyFinalBam.copyFinalBam as copyFinalBam {
            input :
            stepName = stepName + '-copyBam',
            inputBam = bamclipper.bamFile,
            outputBam = finalOutput,
            logFile = logFile
        }
    }

    # Outputs that will be retained when execution is complete
    output {
        String finalFile = select_first([copyFinalBam.copyOutput, pseudoFinal.pseudoFile])
        String bqsrFile = select_first([applyBqsr.bamFile, pseudoBqsr.pseudoFile])
    }
}