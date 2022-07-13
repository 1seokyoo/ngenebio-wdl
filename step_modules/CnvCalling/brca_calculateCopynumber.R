#!/usr/bin/env Rscript

MakeFileSet <- function(path_num){
    sample_sheet <- readLines(paste(path_num,"/SampleList",sep=""),warn=F)
    check_line = 0
    files <- NULL
    for ( i in sample_sheet) {
        files <- append(files,paste(strsplit(i,"\t")[[1]][1],".amplicon.depth.txt",sep=""))
    }

    files
}

MakeCountSet <- function(path_num,files,type,amp_bed){
    tube_list <- TubeNumber(amp_bed)
    amp_count <- length(tube_list[[3]][,1])
    path <- path_num
    total_set <- NULL
    sample_name <- NULL
    sample <- read.delim(paste(path,files[1],sep=""),sep="\t",header=T,stringsAsFactors=F)
    sample_name <- append(sample_name,strsplit(files[1],".amplicon.depth")[[1]][1])
    total_set <- data.frame(round(sample[c(1:amp_count),type],0))
    rownames(total_set) <- sample[c(1:amp_count),1]
    if ( length(files) > 1 ){
        for ( i in 2:length(files) ){
            sample <- read.delim(paste(path,files[i],sep=""),sep="\t",header=T,stringsAsFactors=F)
            sample_name <- append(sample_name,strsplit(files[i],".amplicon.depth")[[1]][1])
            total_set <- cbind(total_set,round(sample[c(1:amp_count),type],0))
        }
    }
    colnames(total_set) <- sample_name

    total_set
}

MakeRefSet <- function(ref_path, type, amp_bed){
    files <- list.files(ref_path, pattern=".amplicon.depth.txt$")

    tube_list <- TubeNumber(amp_bed)
    amp_count <- length(tube_list[[3]][,1])
    path <- ref_path
    total_set <- NULL
    sample_name <- NULL
    sample <- read.delim(paste(path,files[1],sep=""),sep="\t",header=T,stringsAsFactors=F)
    sample_name <- append(sample_name,strsplit(files[1],".amplicon.depth")[[1]][1])
    total_set <- data.frame(round(sample[c(1:amp_count),type],0))
    rownames(total_set) <- sample[c(1:amp_count),1]
    if ( length(files) > 1 ){
        for ( i in 2:length(files) ){
            sample <- read.delim(paste(path,files[i],sep=""),sep="\t",header=T,stringsAsFactors=F)
            sample_name <- append(sample_name,strsplit(files[i],".amplicon.depth")[[1]][1])
            total_set <- cbind(total_set,round(sample[c(1:amp_count),type],0))
        }
    }
    colnames(total_set) <- sample_name

    total_set    

}

MakeMergeSet <- function(ref_set, total_set){
    merge_set <- cbind(ref_set, total_set)

    merge_set

}

TubeNumber <- function(amp_bed){
    
    tube1.num <- NULL
    tube2.num <- NULL
    brca2_count <- 0
    amp_bed_table <- read.delim( amp_bed, sep = "\t", header = F, stringsAsFactors = F )
    for ( i in 1:length(amp_bed_table[,1]) ){
        if ( amp_bed_table[i,9] == 1 ){
            tube1.num <- append(tube1.num,i)
        } else{
            tube2.num <- append(tube2.num,i)
        }
        if ( amp_bed_table[i,1] == "chr13" ){
            brca2_count <- brca2_count + 1
        }
    }

    list(tube1.num, tube2.num, amp_bed_table, brca2_count)

}


Count2Ratio <- function(total_set,amp_bed){
    
    total_set_n <- total_set

    ### Make tube set

    tube_list <- TubeNumber(amp_bed)
    tube1.num <- tube_list[[1]]
    tube2.num <- tube_list[[2]]
    amp.table <- tube_list[[3]]
    brca2_count <- tube_list[[4]]
    brca1_count <- brca2_count + 1
    amp_count <- length(amp.table[,1])
    ### tube 1
    tb1 <- t(total_set_n[c(tube1.num),])
    tb1_n <- as.matrix(tb1)
    for ( i in 1:length(tb1[1,])){
	    total.q <- tb1[,i][tb1[,i] <= quantile(tb1[,i],0.75) & tb1[,i] >= quantile(tb1[,i],0.25)]
	    total.q.index <- which(tb1[,i] <= quantile(tb1[,i],0.75) & tb1[,i] >= quantile(tb1[,i],0.25))
        total.q.median <- median(total.q)
        total.q.sd <- sd(total.q)
        del.index <- NULL
        
        for ( j in 1:length(total.q) ){
            if ( total.q[j] > ( total.q.median + total.q.sd ) | total.q[j] < ( total.q.median - total.q.sd ) ){
                del.index <- append(del.index, j)
            } 
        }
        if ( length(del.index) > 0 ){
            total.q.index <- total.q.index[-del.index]
            total.q <- total.q[-del.index]
        }

	    for ( j in 1:length(tb1[,1])){
		    tb1_n[j,i] <- round( ( tb1[j,i] / sum(tb1[j,]) ) / ( sum(total.q) / sum(tb1[total.q.index,]) ), 3)
	    }
    }

    ### tube 2
    tb2 <- t(total_set_n[c(tube2.num),])
    tb2_n <- as.matrix(tb2)
    for ( i in 1:length(tb2[1,])){
	    total.q <- tb2[,i][tb2[,i] <= quantile(tb2[,i],0.75) & tb2[,i] >= quantile(tb2[,i],0.25)]
	    total.q.index <- which(tb2[,i] <= quantile(tb2[,i],0.75) & tb2[,i] >= quantile(tb2[,i],0.25))
	    total.q.median <- median(total.q)
        total.q.sd <- sd(total.q)
        del.index <- NULL

        for ( j in 1:length(total.q) ){
            if ( total.q[j] > ( total.q.median + total.q.sd ) | total.q[j] < ( total.q.median - total.q.sd ) ){
                del.index <- append(del.index, j)
            }
        }
        if ( length(del.index) > 0 ){
            total.q.index <- total.q.index[-del.index]
            total.q <- total.q[-del.index]
        }

        for ( j in 1:length(tb2[,1])){
		    tb2_n[j,i] <- round( ( tb2[j,i] / sum(tb2[j,]) ) / ( sum(total.q) / sum(tb2[total.q.index,]) ), 3)
	    }
    }

    ### Merge
    total_set_t <- t(total_set_n)
    tb1.idx.count <- 1
    tb2.idx.count <- 1
    for ( i in 1:length(amp.table[,1]) ){
        if ( amp.table[i,9] == 1 ){
            total_set_t[,i] <- tb1_n[,tb1.idx.count]
            tb1.idx.count <- tb1.idx.count + 1
        } else{
            total_set_t[,i] <- tb2_n[,tb2.idx.count]
            tb2.idx.count <- tb2.idx.count + 1
        }
    }

    brca2 <- total_set_t[,c(1:brca2_count)]

    brca1 <- total_set_t[,c(brca1_count:amp_count)]
    
    total_set_final <- t(total_set_t)

    ex.info <- read.delim(amp_bed,sep="\t",header=F,stringsAsFactors=F)

    colnames(ex.info) <- c("chr","start","end","amplicon","exon")
    
    list(brca2,brca1,total_set_final,ex.info)

}


GetUniformity <- function(path_num, total_set, amp_bed){
    
    path <- paste(path_num,"/data/variant/cnv/",sep="")

    total_set_n <- total_set

    ### Make tube set
    
    tube_list <- TubeNumber(amp_bed)
    tube1.num <- tube_list[[1]]
    tube2.num <- tube_list[[2]]

    tb1 <- total_set[c(tube1.num),]
    tb2 <- total_set[c(tube2.num),]

    tb1_uni_mt <- CalUniformity(tb1)
    tb2_uni_mt <- CalUniformity(tb2)

    write.table(tb1_uni_mt,"Uniformity_Tube1.txt", sep="\t", quote=F, row.names = F)
    write.table(tb2_uni_mt,paste(path, "Uniformity_Tube2.txt", sep=""), sep="\t", quote=F, row.names = F)
}

CalUniformity <- function(dt){

    uni_list <- list()
    for ( i in 1:length(colnames(dt)) ) uni_list[[i]] <- c(0,0,0)
    names(uni_list) <- colnames(dt)

    for ( i in ( 1 : length(dt[1,])) ){
	    avg <- mean(dt[,i])
	    avg_05 <- avg*0.5
	    avg_02 <- avg*0.2
	    uni_02_c <- length( which(dt[,i] > avg_02) )
	    uni_05_c <- length( which(dt[,i] > avg_05) )
	    uni_1_c <- length( which(dt[,i] > avg) )
	    total_amp <- length(dt[,1])
	    uni_02_p <- round(uni_02_c / total_amp,2)
	    uni_05_p <- round(uni_05_c / total_amp,2)
	    uni_1_p <- round(uni_1_c / total_amp,2)
	    uni_list[[i]] <- c(uni_1_p,uni_05_p,uni_02_p)
    }   
    
    uni_mt <- do.call(rbind,uni_list)
    uni_out <- cbind(rownames(uni_mt), uni_mt[,c(1:3)])
    colnames(uni_out) <- c("Sample","Uni_1","Uni_0.5","Uni_0.2")
    uni_out
}

TubeThrou <- function(path_num, total_set_raw, ex.info, amp_bed){

    path <- paste(path_num,"/data/variant/cnv/",sep="")

    ### Make tube set
    tube_list <- TubeNumber(amp_bed)
    tube1.num <- tube_list[[1]]
    tube2.num <- tube_list[[2]]
    amp_table <- tube_list[[3]]
    brca2_count <- tube_list[[4]]
    brca1_count <- brca2_count + 1
    amp_count <- length(amp_table[,1])
    ### make amp table
    
    amp_br1 <- apply(ex.info[c(brca1_count:amp_count),],2,rev)
    amp.tmp <- as.matrix(rbind(ex.info[c(1:brca2_count),],amp_br1))
    
    amp <- cbind(amp.tmp[,c(1:4)])
    amp_len <- as.numeric(amp[,3]) - as.numeric(amp[,2])
    
    ### make throuput table

    total_thr <- total_set_raw

    for ( i in (1:length(total_set_raw[,1])) ){
        for ( j in (1:length(total_set_raw[1,])) ){
            total_thr[i,j] <- total_set_raw[i,j] * amp_len[i]
        }
    }

    total_thr_t1 <- total_thr[tube1.num,]
    total_thr_t2 <- total_thr[tube2.num,]
    
    thr_list <- list()
    for ( i in (1:length(total_thr_t1[1,])) ){
        tb1.sum <- round(sum(total_thr_t1[,i])/1000000,1)
        tb2.sum <- round(sum(total_thr_t2[,i])/1000000,1)
        thr_list[[i]] <- c(tb1.sum, tb2.sum)
    }

    thr_mt <- do.call(rbind,thr_list)
    colnames(thr_mt) <- c("Tube1","Tube2")
    rownames(thr_mt) <- colnames(total_thr_t1)
    
    write.table(thr_mt,paste(path,"Throughput_tube.txt", sep=""), sep="\t", quote=F)
}

MakeTable <- function( path_num, brca1, brca2, total_set, total_set_final, ex.info, amp_bed, run_type, outputFile){
    
    sample_name <- colnames(total_set)
    amp_count <- length(total_set[,1])

    tube_list <- TubeNumber(amp_bed)
    tube1.num <- tube_list[[1]]
    tube2.num <- tube_list[[2]]
    amp_table <- tube_list[[3]]
    brca2_count <- tube_list[[4]]
    brca1_count <- brca2_count + 1

    ## raw coverage matrix    
    total_set_br1 <- apply(total_set[c(brca1_count:amp_count),],2,rev)
    total_set_sorted <- as.matrix(rbind(total_set[c(1:brca2_count),],total_set_br1))
    ##
      
    brca1_ratio_high <- NULL
    brca1_ratio_low <- NULL
    for ( i in (1:length(brca1[1,]))){
            c.pos <- i
            tmp_set <- brca1[,c.pos][brca1[,c.pos] <= quantile(brca1[,c.pos],0.75) & brca1[,c.pos] >= quantile(brca1[,c.pos],0.25)]
            tmp_set.index <- which(brca1[,c.pos] <= quantile(brca1[,c.pos],0.75) & brca1[,c.pos] >= quantile(brca1[,c.pos],0.25))
            
            tmp_list <- outFilter(tmp_set, tmp_set.index)
            tmp_set <- tmp_list[1][[1]]

            brca1_ratio_high <- append(brca1_ratio_high, round(max(tmp_set)*1.1,3))
            brca1_ratio_low <- append(brca1_ratio_low, round(min(tmp_set)*0.9,3))
    }


    brca2_ratio_high <- NULL
    brca2_ratio_low <- NULL
    for ( i in (1:length(brca2[1,]))){
            c.pos <- i
            tmp_set <- brca2[,c.pos][brca2[,c.pos] <= quantile(brca2[,c.pos],0.75) & brca2[,c.pos] >= quantile(brca2[,c.pos],0.25)]
            tmp_set.index <- which(brca2[,c.pos] <= quantile(brca2[,c.pos],0.75) & brca2[,c.pos] >= quantile(brca2[,c.pos],0.25))

            tmp_list <- outFilter(tmp_set, tmp_set.index)
            tmp_set <- tmp_list[1][[1]]

            brca2_ratio_high <- append(brca2_ratio_high, round(max(tmp_set)*1.1,3))
            brca2_ratio_low <- append(brca2_ratio_low, round(min(tmp_set)*0.9,3))
    }

    total_high <- c(brca2_ratio_high,rev(brca1_ratio_high))
    total_low <- c(brca2_ratio_low,rev(brca1_ratio_low))

    ## data list set
    tmp.dt.final.list <- list()
    
    inter_list <- list()
    
    raw_dp_list <- list()

    
    ## Cal dist val

    for ( i in (1:length(total_set_final[,1])) ){
        c.pos <- i
        tmp.q <- total_set_final[c.pos,][total_set_final[c.pos,] <= quantile(total_set_final[c.pos,],0.75) & total_set_final[c.pos,] >= quantile(total_set_final[c.pos,],0.25)]
        tmp.q.index <- which(total_set_final[c.pos,] <= quantile(total_set_final[c.pos,],0.75) & total_set_final[c.pos,] >= quantile(total_set_final[c.pos,],0.25))
        sample.num <- length(tmp.q)
        sample.mean <- mean(tmp.q)
        sample.sd <- sd(tmp.q)
        error <- qnorm(0.999)*sample.sd*2/sqrt(sample.num)
        inter.left <- sample.mean - error
        inter.right <- sample.mean + error

        inter.low <- round(qnorm(0.001, mean = sample.mean, sd = sample.sd*2 ), 3)
        inter.high <- round(qnorm(0.999, mean = sample.mean, sd = sample.sd*2 ), 3)

        raw_dp.median <- round(median(total_set_sorted[c.pos,tmp.q.index]))
        raw_dp.mean <- round(mean(total_set_sorted[c.pos,tmp.q.index]))

        sample.val.total <- NULL
        sample.val.total.raw <- NULL
        for ( t in (1:length(total_set_sorted[1,])) ){
            sample.val <- total_set_final[c.pos,t]
            sample.val.raw <- total_set_sorted[c.pos,t]
            sample.val.total <- append(sample.val.total, sample.val)
            sample.val.total.raw <- append(sample.val.total.raw, sample.val.raw)
        }

        inter_list[[i]] <- c(inter.low, inter.high, sample.val.total)
        raw_dp_list[[i]] <- c(raw_dp.median, raw_dp.mean, sample.val.total.raw)
    }

    inter_mt <- do.call(rbind,inter_list)
    raw_dp.mt <- do.call(rbind,raw_dp_list)

    ## Write table

    for ( k in 1:length(sample_name) ){
        if ( run_type > 1 & k != length(sample_name) ) next

        tmp.dt.final <- cbind(ex.info[,c(6,5,7)], total_low, total_high, total_set_final[,k])

        colnames(tmp.dt.final) <- c("Gene","Exon","Amplicon","Min","Max","Ratio")
        write.table(tmp.dt.final,paste(path_num,sample_name[k],".copynumber.txt",sep=""),sep="\t",quote=F,row.names = F)
        writeLines(c("Complete to ", sample_name[k]), outputFile)
    }

}

outFilter <- function(set,set.index){

    set.median <- median(set)
    set.sd <- sd(set)
    del.index <- NULL

    for ( i in 1:length(set) ){
        if ( set[i] > ( set.median + set.sd ) | set[i] < ( set.median - set.sd ) ){
            del.index <- append(del.index, i)
        }   
    }

    if ( length(del.index) > 0 ){
        set.index <- set.index[-del.index]
        set <- set[-del.index]
    }

    list(set, set.index)

}

RefSetCNV <- function(ref_path, total_set, amp_bed, path_num, outputFile){
    ref_set <- MakeRefSet(ref_path, 3, amp_bed)
    for ( i in 1:length(total_set[1,]) ){
        cnv_set <- MakeMergeSet(ref_set, total_set[,i])
        old_name <- colnames(cnv_set)
        old_name[length(ref_set[1,])+1] <- colnames(total_set[i])
        colnames(cnv_set) <- old_name
        brca_list <- Count2Ratio(cnv_set,amp_bed)
        brca2 <- brca_list[1][[1]]
        brca1 <- brca_list[2][[1]]
        total_set_final <- brca_list[3][[1]]
        ex.info <- brca_list[4][[1]]
        MakeTable(path_num,brca1,brca2,cnv_set,total_set_final,ex.info,amp_bed,2,outputFile)
    }

}

RunSetCNV <- function(total_set, amp_bed, path_num, outputFile){
    brca_list <- Count2Ratio(total_set,amp_bed)
    brca2 <- brca_list[1][[1]]
    brca1 <- brca_list[2][[1]]
    total_set_final <- brca_list[3][[1]]
    ex.info <- brca_list[4][[1]]
    MakeTable(path_num,brca1,brca2,total_set,total_set_final,ex.info,amp_bed,1,outputFile)
}

RunCNV <- function(path_num, amp_bed, ref_path, total_set, outputFile){
    if ( length(total_set[1,]) < 6 ){
        RefSetCNV(ref_path, total_set, amp_bed, path_num, outputFile)
    } else{
        RunSetCNV(total_set, amp_bed, path_num, outputFile)
    }
}


args <- (commandArgs(TRUE))
path_num <- args[1]
ref_path <- args[2]
amp_bed <- args[3]
outputFile <- file(args[4])


# files <- MakeFileSet(path_num)
files <- list.files(path_num, pattern=".amplicon.depth.txt$")

total_set <- MakeCountSet(path_num,files,3,amp_bed)

RunCNV(path_num, amp_bed, ref_path, total_set, outputFile)

close(outputFile)