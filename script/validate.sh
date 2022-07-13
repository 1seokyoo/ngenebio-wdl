#!/bin/bash
#================================================================================
#
# FILE: validate.sh
#
# DESCRIPTION: System validation
# Note:
# REQUIREMENTS:
# BUGS:
# NOTES:
# AUTHOR: Chnagbum Hong cbhong@ngenebio.com
# COMPANY: ngenebio
# VERSION: 1.0
# CREATED: 2015_12_1
# REVISION:
#================================================================================

script_name()
{
   script_name=`echo $0 | sed 's/-//'`
   echo `basename $script_name`
}

timestamp()
{
   date +"%T %D"
}

#General logging function not to be called on itw own
log()
{
   message=$1
   # echo -e "[$(timestamp)]" $inputStep "["$message >> $logFile 
   echo -e $inputStep "[$(timestamp)]" "["$message >> $logFile 
}

log_checkpoint()
{
   log "CHECKPOINT] > $1"
}

log_progress()
{
   # log "PROGRESS] > $1"
   echo -e "[PROGRESS]" $inputStep "[$(timestamp)] > " $1 >> $logFile 

}

log_task()
{
   # log "TASK] > $1"
   echo -e "[TASK]" $inputStep "[$(timestamp)] > " $1 >> $logFile 

}

log_error()
{
   # log "ERROR] > $1"
   echo -e "[ERROR]" $inputStep "[$(timestamp)] > " $1 >> $logFile 
   exit 1
}

log_info()
{
   # log "INFO] > $1"
   echo -e "[INFO]" $inputStep "[$(timestamp)] > " $1 >> $logFile 

}

log_cmd()
{
   # log "CMD] > $1"
   echo -e "[CMD]" $inputStep "[$(timestamp)] > " $1 >> $logFile 
}

#Wrap 3rd party tool calls so we can
#log how they were executed
wrap()
{
	log_cmd "$*" 
	eval $* 2>> $logFile
}

md5sum_check()
{
    if [[ -f ${1} && -f "${1}.md5" && `md5sum ${1} | cut -f1 -d" "` = `cat "${1}.md5" | cut -f1 -d" "` ]]; then
        echo 'True'
    else
        echo 'False'
    fi
}

vcf_check()
{
   if [[ -f $1 && `echo $1 | grep "\.gz"` = $1 ]]; then
      local variant_count=`zcat $1 | grep -v "^#" | wc -l`
   elif [[ -f $1 && `echo $1 | grep "\.gz"` != $1 ]]; then
      local variant_count=`cat $1 | grep -v "^#" | wc -l`
   else
      local variant_count=0
   fi

   if [[ ${variant_count} > 0 ]]; then
      echo 'True'
   else
      echo 'False'
   fi
}