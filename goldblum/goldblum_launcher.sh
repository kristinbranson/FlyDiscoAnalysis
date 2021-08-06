#! /bin/bash

escaped_bash_profile_path=${1}
escaped_fly_disco_analysis_folder_path=${2}
pi_last_name=${3}
escaped_goldblum_logs_folder_path=${4}
date_as_string=`date +%Y-%m-%d`
goldblum_log_file_name="goldblum-${date_as_string}.log"
goldblum_log_file_path="${escaped_goldblum_logs_folder_path}/${goldblum_log_file_name}"

. /misc/lsf/conf/profile.lsf
. ${escaped_bash_profile_path}
cd ${escaped_fly_disco_analysis_folder_path}
bsub -n1 -P ${pi_last_name} -o ${goldblum_log_file_path} -e ${goldblum_log_file_path} \
  /misc/local/matlab-2019a/bin/matlab -nodisplay -batch 'modpath; goldblum(true, true);'
