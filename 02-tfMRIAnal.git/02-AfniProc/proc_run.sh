#!/bin/bash
# Copyright (c) 2022 Institute of Health and Medical Technology, Hefei Institutes of Physical Science, CAS, All Rights Reserved.
# @Time     : 2022/12/20 13:09:12
# @Author   : ZL.Z
# @Email    : zzl1124@mail.ustc.edu.cn
# @Reference: None
# @FileName : proc_run.sh
# @Software : Shell; VS Code; Ubuntu 18.04.5 LTS (GNU/Linux 5.4.0-59-generic x86_64)
# @Hardware : Intel Xeon 6258R 2.7G 28C*2
# @Version  : V1.0.0 - ZL.Z: 2022/12/20
#             First version.
# @License  : None
# @Brief    : 执行AFNI的预处理及GLM：afni_proc.py

beginTime=$(date +%s)
# get basement path of the corrent shell file
# sh_basePath=$(cd `dirname $0`;pwd)
sh_basePath="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# 将中文路径设置软链接
# 判断主文件夹下是否存在上次设置的软链接，存在则删除
if [[ -L ~/temp_ln_src_proc ]]; then
	rm -rf ~/temp_ln_src_proc
fi
ln -s "/home/medicaldata2/ZZLData/Datasets/audio/言语与认知老化实验/言语与认知老化" ~/temp_ln_src_proc # 创建软链接
# source data path
src_dat_dir=~/temp_ln_src_proc/data/raw_data/t-fMRI/MRI/NIfTI/
# output path
subj_dirs=~/temp_ln_src_proc/data/raw_data/t-fMRI/MRI/afni_proc_res/
# timing files path
timing_dirs=~/data/言语与认知老化实验/言语与认知老化/analysis/MRI_process/func/CreatingTimingFilesAFNI/timing_file/dmBLOCK/

subid_for_run="20230227087 20230308093 20230213064 20230424102 20230222079 20230220076"

# parallel running
for subj_dir in "$src_dat_dir"/*; do
	{
		subid=$(basename "$subj_dir")
		# if [[ $subid_for_run =~ $subid ]]; then
		if [[ $subid != -all* && $subid != -func* && ! -d $subj_dirs/"$subid".results && -d $timing_dirs/"$subid" ]]; then
			echo '----------------------------' "$subid" '----------------------------'
			tcsh "$sh_basePath"/proc.default_subid "$subj_dir" $timing_dirs/"$subid" $subj_dirs/"$subid".results
		fi
	} &
done

# Wait for all background tasks to complete before proceeding to the next step
wait
if [[ -L ~/temp_ln_src_proc ]]; then
	rm -rf ~/temp_ln_src_proc # 删除创建的软链接
fi
endTime=$(date +%s)
echo "Finished in" $((endTime - beginTime)) "s"
