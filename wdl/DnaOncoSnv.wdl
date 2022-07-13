#================================================================================
#
# FILE: DnaOncoSnv.wdl
#
# DESCRIPTION: NGeneBio DNA ONCOaccuPanel DnaOncoSnv workflow
# Note:
# REQUIREMENTS:
# BUGS:
# NOTES:
# AUTHOR: Juyeon hong hong.juyeon@ngenebio.com
# COMPANY: ngenebio
# VERSION: 0.1
# CREATED: 2022.05.19
# REVISION:
# RUN COMMAND: /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/cromwell-50.jar run /opt/ngenebio/pipeline/wdl/DnaOncoSnv.wdl -i /opt/ngenebio/test/testOncoSnv.json
# VALIDATION : /opt/ngenebio/app/java/jre1.8.0_144/bin/java -Dconfig.file=/opt/ngenebio/pipeline/cromwell.conf -jar /opt/ngenebio/app/cromwell/womtool-50.jar validate /opt/ngenebio/pipeline/wdl/DnaOncoSnv.wdl -i /opt/ngenebio/test/testOncoSnv.json
#================================================================================

import '/opt/ngenebio/pipeline/step_modules/SnvCalling/VardictV16Run.wdl' as Vardict
import '/opt/ngenebio/pipeline/step_modules/SnvCalling/VcfFilter.wdl' as VcfFilter
import '/opt/ngenebio/pipeline/step_modules/SnvCalling/VcfValidate.wdl' as VcfValidation
import '/opt/ngenebio/pipeline/step_modules/SnvCalling/RunMutect2.wdl' as Mutect2
import '/opt/ngenebio/pipeline/step_modules/SnvCalling/WhitelistProcess.wdl' as WhitelistProcess
import '/opt/ngenebio/pipeline/step_modules/SnvCalling/VcfMerge.wdl' as VcfMerge
import '/opt/ngenebio/pipeline/step_modules/SnvCalling/Vep.wdl' as Vep
import '/opt/ngenebio/pipeline/wdl/Setting.wdl' as Setting

workflow DnaOncoSnv {
    
    # Input Parameter
    String stepName = 'SnvIndelCalling'
    String logDir
    String resultDir
    String variantDir 
    String assayDir 
    String alignDir
    String tmpDir

    String inputFinalBam 
    String sampleName
    Int sampleId
    Int runId
    String sampleDisease
    
    String inputFasta
    String inputTargetAllGene  
    String inputTargetHrd
    String inputTargetWhitelist
    String inputTargetWhitelistTxt 
    String inputSnvFormat
    Int cpu = 2
    Int ram = 4
    
    String logFile = logDir + sampleId + '.' + runId + '.log.snv.stderr.txt'
    String finalVardictAllGeneOutput = variantDir + sampleName + '.'+ 'all' + '.final.vcf'
    String finalVardictHrdOutput = variantDir + sampleName + '.'+ 'hrd' + '.final.vcf'
    String finalVardictWhitelistOutput = variantDir + sampleName + '.'+ 'whitelist' + '.final.vcf'
    String finalVardictWhitelistInfo = variantDir + sampleName + '.'+ 'whitelist' + '.variant.txt'
    String finalMutect2AllGeneOutput = variantDir + sampleName + '.'+ 'mutect' + '.final.vcf'
    String finalMergedOutput = variantDir + sampleName + '.total.raw.vcf'
    String finalVcfOutput = resultDir + sampleName + '.snv.vcf'


    ##### Variant All Gene #####
    String stepNameAll = stepName + '-All'
    call Setting.checkMd5Sum as checkVardictAllGene {
        input:
            stepName = stepNameAll + '-checkMd5Sum',
            logFile = logFile,
            inputFile = finalVardictAllGeneOutput
    }

    if (checkVardictAllGene.workMd5Sum == 'true') {
        call Setting.pseudoOutput as vardictAllGenePseudoOutput {
            input:
                inputFile = finalVardictAllGeneOutput
        }
    }
    
    if (checkVardictAllGene.workMd5Sum != 'true') {

        call Vardict.vardictCall as runVardictAllGene {
            input:
                stepName = stepNameAll + '-vardictCall',
                inputTarget = inputTargetAllGene,
                inputBam = inputFinalBam, 
                sampleName = sampleName,
                inputFasta = inputFasta,
                mode = 'all',
                outputDir = variantDir,
                optionFreq = 0.01,
                optionGenomicCall = false,
                logFile = logFile,
                cpu = cpu
        }

        call Vardict.teststrandbias as runStrandbiasAllGene {
            input:
                stepName = stepNameAll + '-teststrandbias',
                inputVardict = runVardictAllGene.rawVardictFile,
                outputDir = variantDir,
                sampleName = sampleName,
                mode = 'all',
                logFile = logFile,
                cpu = cpu
        }

        call Vardict.var2Vcf as runVar2VcfAllGene {
            input:
                stepName = stepNameAll + '-var2vcf',
                inputStrandbias = runStrandbiasAllGene.rawStrandbiasFile,
                optionAo = 3,
                optionFreq = 0.03,
                outputDir = variantDir,
                sampleName = sampleName,
                mode = 'all',
                logFile = logFile,
                cpu = cpu
        }

        call VcfValidation.vcfValidate as validateVcfAllGene {
            input:
                stepName = stepNameAll + '-validation',
                inputVcf = runVar2VcfAllGene.rawVcfFile,
                inputFasta = inputFasta,
                sampleName = sampleName,
                outputDir = variantDir,
                mode = 'all',
                tmpDir = tmpDir,
                logFile = logFile,
                cpu = cpu,
                ram = ram
        }

        call VcfFilter.vcfFilter as filterVcfAllGene {
            input:
                stepName = stepNameAll + '-vcfFilter',
                inputVcf = validateVcfAllGene.formatVcfFile,
                sampleName = sampleName,
                outputVcf = finalVardictAllGeneOutput,
                databaseList = 'dbSNPBuildID,gnomAD_AF,PON_DB,KRG_DB,BLACKLIST_DB,ExAC_EAS_AF',
                vcfFormat = 'DP',
                cutoff = 30,
                clinvarFilter = 'None',
                logFile = logFile
        }
    }

    String vardictAllVcfFile = select_first([filterVcfAllGene.finalVcfFile, vardictAllGenePseudoOutput.pseudoFile])


    ##### VarDict HRD Gene #####
    String stepNameHrd = stepName + '-HRD'
    call Setting.checkMd5Sum as checkVardictHrd {
        input:
            stepName = stepNameHrd + '-checkMd5Sum',
            logFile = logFile,
            inputFile = finalVardictHrdOutput
    }

    if (checkVardictHrd.workMd5Sum == 'true') {
        call Setting.pseudoOutput as vardictHrdPseudoOutput{
            input:
                inputFile = finalVardictHrdOutput
        }
    }
    
    if (checkVardictHrd.workMd5Sum != 'true') {

        call Vardict.vardictCall as runVardictHrd {
            input:
                stepName = stepNameHrd + '-vardictCall',
                inputTarget = inputTargetHrd,
                inputBam = inputFinalBam, 
                sampleName = sampleName,
                inputFasta = inputFasta,
                mode = 'hrd',
                outputDir = variantDir,
                optionFreq = 0.1,
                optionGenomicCall = false,
                logFile = logFile,
                cpu = cpu
        }

        call Vardict.teststrandbias as runStrandbiasHrd {
            input:
                stepName = stepNameHrd + '-teststrandbias',
                inputVardict = runVardictHrd.rawVardictFile,
                outputDir = variantDir,
                sampleName = sampleName,
                mode = 'hrd',
                logFile = logFile,
                cpu = cpu
        }

        call Vardict.var2Vcf as runVar2VcfHrd {
            input:
                stepName = stepNameHrd + '-var2vcf',
                inputStrandbias = runStrandbiasHrd.rawStrandbiasFile,
                optionAo = 10,
                optionFreq = 0.1,
                outputDir = variantDir,
                sampleName = sampleName,
                mode = 'hrd',
                logFile = logFile,
                cpu = cpu
        }

        call VcfValidation.vcfValidate as validateVcfHrd {
            input:
                stepName = stepNameHrd + '-validation',
                inputVcf = runVar2VcfHrd.rawVcfFile,
                inputFasta = inputFasta,
                sampleName = sampleName,
                outputDir = variantDir,
                mode = 'hrd',
                tmpDir = tmpDir,
                logFile = logFile,
                cpu = cpu,
                ram = ram
        }

        call VcfFilter.vcfFilter as filterVcfHrd {
            input:
                stepName = stepNameHrd + '-vcfFilter',
                inputVcf = validateVcfHrd.formatVcfFile,
                sampleName = sampleName,
                outputVcf = finalVardictHrdOutput,
                databaseList = 'BLACKLIST_DB',
                vcfFormat = 'None',
                cutoff = 'None',
                clinvarFilter = 'Yes',
                logFile = logFile
        }
    }

    String vardictHrdVcfFile = select_first([filterVcfHrd.finalVcfFile,
                                             vardictHrdPseudoOutput.pseudoFile])


    ##### Variant Whitelist Gene #####
    String stepNameWhi = stepName + '-WHI'
    call Setting.checkMd5Sum as checkVardictWhitelist {
        input:
            stepName = stepNameWhi + '-checkMd5Sum',
            logFile = logFile,
            inputFile = finalVardictWhitelistOutput
    }

    if (checkVardictWhitelist.workMd5Sum == 'true') {
        call Setting.pseudoOutput as vardictWhitelistPseudoOutput{
            input:
                inputFile = finalVardictWhitelistOutput
        }
        call Setting.pseudoOutput as vardictWhitelistPseudoInfo{
            input:
                inputFile = finalVardictWhitelistInfo
        }
    }

    if (checkVardictWhitelist.workMd5Sum != 'true') {

        call WhitelistProcess.makeTarget as makeWhitelistTarget {
            input:
                stepName = stepNameWhi + '-makeTarget',
                sampleName = sampleName,
                logFile = logFile,
                inputDisease = sampleDisease,
                inputTarget = inputTargetWhitelistTxt,
                outputDir = assayDir,
        }

        call WhitelistProcess.checkTarget as checkWhitelistTarget{
            input:
                stepName = stepNameWhi + '-checkTarget',
                inputTarget = makeWhitelistTarget.listFile,
                outputVcf = finalVardictWhitelistOutput,
                logFile = logFile
        }
        
        if (checkWhitelistTarget.finalTxt != 'true' ) {
            call Setting.pseudoOutput as vardictNonWhitelistPseudoOutput{
                input:
                    inputFile = checkWhitelistTarget.finalVcfFile
            }

            call Setting.pseudoOutput as vardictNonWhitelistPseudoInfo{
                input:
                    inputFile = checkWhitelistTarget.finalVariantFile
            }
        }

        if (checkWhitelistTarget.finalTxt == 'true' ) {

            call Vardict.vardictCall as runVardictWhitelist {
                input:
                    stepName = stepNameWhi+ '-vardictCall',
                    inputTarget = inputTargetWhitelist,
                    inputBam = inputFinalBam, 
                    sampleName = sampleName,
                    inputFasta = inputFasta,
                    mode = 'whitelist',
                    outputDir = variantDir,
                    optionFreq = 0.001,
                    optionGenomicCall = true,
                    logFile = logFile,
                    cpu = cpu
            }

            call Vardict.teststrandbias as runStrandbiasWhitelist {
                input:
                    stepName = stepName + '-WHI' + '-teststrandbias',
                    inputVardict = runVardictWhitelist.rawVardictFile,
                    outputDir = variantDir,
                    sampleName = sampleName,
                    mode = 'whitelist',
                    logFile = logFile,
                    cpu = cpu
            }

            call Vardict.var2Vcf as runVar2VcfWhitelist {
                input:
                    stepName = stepName + '-WHI' + '-var2vcf',
                    inputStrandbias = runStrandbiasWhitelist.rawStrandbiasFile,
                    optionAo = 3,
                    optionFreq = 0.001,
                    outputDir = variantDir,
                    sampleName = sampleName,
                    mode = 'whitelist',
                    logFile = logFile,
                    cpu = cpu
            }

            call WhitelistProcess.whitelistFilter as filterWhitelist {
                input:
                    stepName = stepName + '-WHI' + '-whitelistFilter',
                    logFile = logFile,
                    inputVcf = runVar2VcfWhitelist.rawVcfFile,
                    inputTarget = makeWhitelistTarget.listFile
            }

            call VcfValidation.vcfValidate as validateVcfWhitelist {
                input:
                    stepName = stepName + '-WHI' + '-validation',
                    inputVcf = filterWhitelist.vcfFile,
                    inputFasta = inputFasta,
                    sampleName = sampleName,
                    outputDir = variantDir,
                    mode = 'whitelist',
                    tmpDir = tmpDir,
                    logFile = logFile,
                    cpu = cpu,
                    ram = ram
            }

            call WhitelistProcess.filterDepth as filterDepthVcfWhitelist {
                input:
                    stepName = stepName + '-WHI' + '-filterDepth',
                    logFile = logFile,
                    inputVcf = validateVcfWhitelist.formatVcfFile,
                    vcfFormat = 'VD',
                    cutoff = 3
            }

    
            call WhitelistProcess.whitelistInfo as runInfoWhitelist {
                input:
                    stepName = stepName + '-WHI' + '-whitelistinfo',
                    inputVcf = filterDepthVcfWhitelist.vcfFile,
                    logFile = logFile,
                    inputTxt = filterWhitelist.txtFile
            }
            
            # call WhitelistProcess.vep as runVepVcfWhitelist {
            #     input:
            #         inputVcf = filterDepthVcfWhitelist.vcfFile,
            #         inputFasta = inputFasta,
            #         outputVcf = finalVardictWhitelistOutput,
            #         logFile = logFile
            # }

        }
    }
    
    String vardictWhiVcfFile = select_first([filterDepthVcfWhitelist.vcfFile,
                                             vardictNonWhitelistPseudoOutput.pseudoFile,
                                             vardictWhitelistPseudoOutput.pseudoFile])


    ##### Mutect All Gene #####
    String stepNameMut = stepName + '-MUT'
    call Setting.checkMd5Sum as checkMutect2AllGene {
        input:
            stepName = stepNameMut + '-checkMd5Sum',
            logFile = logFile,
            inputFile = finalMutect2AllGeneOutput,
    }

    if (checkMutect2AllGene.workMd5Sum == 'true') {
        call Setting.pseudoOutput as mutect2AllGenePseudoOutput{
            input:
                inputFile = finalMutect2AllGeneOutput
        }
    }

    if (checkMutect2AllGene.workMd5Sum != 'true') {

        call Mutect2.cleanSam as runCleanSamMutect2 {
            input:
                stepName = stepNameMut + '-cleanSam',
                inputBam = inputFinalBam,
                sampleName = sampleName,
                outputDir = alignDir,
                logFile = logFile,
                tmpDir = tmpDir,
                cpu = cpu,
                ram = ram
        }

        call Mutect2.mutect2 as callMutect2 {
            input:
                stepName = stepNameMut + '-Mutect2',
                inputBam = runCleanSamMutect2.bamFile,
                inputTarget = inputTargetAllGene ,
                inputFasta = inputFasta,
                sampleName = sampleName,
                outputDir = variantDir,
                logFile = logFile,
                tmpDir = tmpDir,
                cpu = cpu,
                ram = ram
        }


        call Mutect2.mutect2PostProcessing as postProcessMutect2 {
            input:
                stepName = stepNameMut + '-PostProcessing',
                inputBam = runCleanSamMutect2.bamFile,
                inputVcf = callMutect2.vcfFile,
                inputTarget = inputTargetAllGene,
                inputFasta = inputFasta,
                inputF1R2 = callMutect2.f1R2File,
                sampleName = sampleName,
                outputDir = variantDir,
                tmpDir = tmpDir,
                logFile = logFile,
                cpu = cpu,
                ram = ram
        }

        
        call VcfValidation.vcfValidate as validateVcfMutect2 {
            input:
                stepName = stepNameMut + '-validation',
                inputVcf = postProcessMutect2.vcfFile,
                inputFasta = inputFasta,
                sampleName = sampleName,
                outputDir = variantDir,
                mode = 'mutect2',
                tmpDir = tmpDir,
                logFile = logFile,
                cpu = cpu,
                ram = ram
        }

        call VcfFilter.vcfFilter as filterVcfMutect2 {
            input:
                stepName = stepNameMut + '-vcfFilter',
                inputVcf = validateVcfMutect2.formatVcfFile,
                sampleName = sampleName,
                outputVcf = finalMutect2AllGeneOutput,
                databaseList = 'dbSNPBuildID,gnomAD_AF,PON_DB,KRG_DB,BLACKLIST_DB,ExAC_EAS_AF',
                vcfFormat = 'None',
                cutoff = 'None',
                clinvarFilter = 'None',
                logFile = logFile
        }
    }

    String mutectVcfFile = select_first([filterVcfMutect2.finalVcfFile, mutect2AllGenePseudoOutput.pseudoFile])


    ##### VCF Merge #####
    String stepNameMerge = stepName + '-Merge'
    call Setting.checkMd5Sum as checkMergedVcf {
        input:
            stepName = stepNameMerge + '-checkMd5Sum',
            logFile = logFile,
            inputFile = finalMergedOutput,
    }

    if (checkMergedVcf.workMd5Sum == 'true') {
        call Setting.pseudoOutput as mergeVcfPseudoOutput{
            input:
                inputFile = finalMergedOutput
        }
    }

    if (checkMergedVcf.workMd5Sum != 'true') {
        call VcfMerge.vcfMerge as mergeVcf {
            input:
                stepName = stepNameMerge + '-vcfMerge',
                inputAllGeneVcf = vardictAllVcfFile,
                inputHrdVcf = vardictHrdVcfFile,
                inputWhitelistVcf = vardictWhiVcfFile,
                inputMutectVcf =  mutectVcfFile,
                inputSnvFormat = inputSnvFormat,
                sampleName = sampleName,
                outputDir = variantDir,
                outputVcf = finalMergedOutput,
                logFile = logFile
        } 
    }

    String mergedVcfFile = select_first([mergeVcf.vcfFile, mergeVcfPseudoOutput.pseudoFile])


    ##### VEP Filter #####
    String stepNameVep = stepName + '-VepFilter'
    call Setting.checkMd5Sum as checkVepVcf {
        input:
            stepName = stepNameVep + '-checkMd5Sum',
            logFile =logFile,
            inputFile = finalVcfOutput
    }

    if (checkVepVcf.workMd5Sum == 'true') {
        call Setting.pseudoOutput as vepVcfPseudoOutput{
            input:
                inputFile = finalVcfOutput
        }
    }

    if (checkVepVcf.workMd5Sum != 'true') {
        call Vep.run as runVep {
            input :
                stepName = stepNameVep + '-run',
                inputVcf = mergedVcfFile,
                inputFasta = inputFasta,
                sampleName = sampleName,
                outputDir = variantDir,
                ensemblVersion = 92,
                logFile = logFile
        }

        call Vep.filter as runVepFilter {
            input :
                stepName = stepNameVep + '-filter',
                inputVcf = runVep.vcfFile,
                sampleName = sampleName,
                outputDir = variantDir,
                logFile = logFile                
        }

        call Setting.copyFinalFile {
            input :
                stepName = stepName + '-copyFinalFile',
                logFile = logFile,
                inputFile = runVepFilter.vcfFile,
                outputName = finalVcfOutput
        }
    }

    output {
        String finalVcfFile = select_first([copyFinalFile.outputFile, vepVcfPseudoOutput.pseudoFile])
        String whitelistInfoFile = select_first([runInfoWhitelist.txtFile,
                                                 vardictNonWhitelistPseudoInfo.pseudoFile,
                                                 vardictWhitelistPseudoInfo.pseudoFile])
    }
}