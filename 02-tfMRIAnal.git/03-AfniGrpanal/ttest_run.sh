#!/bin/bash
# Copyright (c) 2023 Institute of Health and Medical Technology, Hefei Institutes of Physical Science, CAS, All Rights Reserved.
# @Time     : 2023/8/17
# @Author   : ZL.Z
# @Email    : zzl1124@mail.ustc.edu.cn
# @Reference: None
# @FileName : ttest_run.sh
# @Software : Shell; VS Code; Ubuntu 18.04.5 LTS (GNU/Linux 5.4.0-59-generic x86_64)
# @Hardware : Intel Xeon 6258R 2.7G 28C*2
# @Version  : V1.0.0 - ZL.Z: 2023/8/17
#             First version.
# @License  : None
# @Brief    : 执行AFNI的t-test组分析：uber_ttest.py

beginTime=$(date +%s)
begin_time=$(date)
# get basement path of the corrent shell file
# sh_basePath=$(cd `dirname $0`;pwd)
sh_basePath="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# 将中文路径设置软链接
# 判断主文件夹下是否存在上次设置的软链接，存在则删除
if [[ -L ~/temp_ln_src_ttest ]]; then
    rm -rf ~/temp_ln_src_ttest
fi
ln -s "/home/medicaldata2/ZZLData/Datasets/audio/言语与认知老化实验/言语与认知老化" ~/temp_ln_src_ttest # 创建软链接
# proc result path
proc_res_dirs=~/temp_ln_src_ttest/data/raw_data/t-fMRI/MRI/afni_proc_res/

if [[ -L ~/temp_ln_subj_ttest ]]; then
    rm -rf ~/temp_ln_subj_ttest
fi
ln -s "$sh_basePath" ~/temp_ln_subj_ttest # 创建软链接
# output path
output_dir=~/temp_ln_subj_ttest/results/
# Brainnetome atlas: Broca(29,33,35,37,39);PMC(53,54,57,58,59,60,61,62,67,68);dlPFC(5,6,11,12,15,16,19,20,21,22);
# ACC(177,178,179,180,183,184,187,188);
# hippo(15,216,217,218);thalamus(231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246)
3dcalc -a roi_mask/brainnetome/BN_Atlas_246_1mm.nii.gz \
-expr "amongst(a,29,33,35,37,39,53,54,57,58,59,60,61,62,67,68,5,6,11,12,15,16,19,20,21,22,177,178,179,180,183,184,187,188,
215,216,217,218,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246)" \
       -prefix roi_mask/roi_all.nii.gz -overwrite
3dresample -master $proc_res_dirs/20220928002.results/stats.20220928002_REML+tlrc -input roi_mask/roi_all.nii.gz \
           -prefix roi_mask/roi_all_mask.nii.gz -overwrite
3dcalc -a roi_mask/brainnetome/BN_Atlas_246_1mm.nii.gz \
-expr "amongst(a,29,33,35,37,39,53,54,57,58,59,60,61,62,67,68,215,216,217,218,
231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246)" \
       -prefix roi_mask/roi_broca_pmc_hip_tha.nii.gz -overwrite
3dresample -master $proc_res_dirs/20220928002.results/stats.20220928002_REML+tlrc -input roi_mask/roi_broca_pmc_hip_tha.nii.gz \
           -prefix roi_mask/roi_broca_pmc_hip_tha_mask.nii.gz -overwrite
3dcalc -a roi_mask/brainnetome/BN_Atlas_246_1mm.nii.gz \
-expr "amongst(a,5,6,11,12,15,16,19,20,21,22,177,178,179,180,183,184,187,188)" \
       -prefix roi_mask/roi_dlpfc_acc.nii.gz -overwrite
3dresample -master $proc_res_dirs/20220928002.results/stats.20220928002_REML+tlrc -input roi_mask/roi_dlpfc_acc.nii.gz \
           -prefix roi_mask/roi_dlpfc_acc_mask.nii.gz -overwrite

use_reml="_REML"  # using REMLfit ("_REML") or not ("")
# labels="base dt1 dt2 dt3 dt4 dt1-base dt2-base dt1-dt3 dt2-dt4 dt1-base-dt3 dt2-base-dt4 dt1-dt3-base dt2-dt4-base"
labels="dt1 dt2 dt1-base dt2-base dt1-dt3 dt2-dt4"

# parallel running
for key in $labels; do
    {
        echo '----------------------------' T-test: "$key" '----------------------------'
        mask_f=~/temp_ln_subj_ttest/roi_mask/roi_all_mask.nii.gz
        if [[ $key == "dt1" || $key == "dt2" ]]; then
            mask_f=~/temp_ln_subj_ttest/roi_mask/roi_broca_pmc_hip_tha_mask.nii.gz
        elif [[ $key == "dt1-base" || $key == "dt2-base" || $key == "dt1-dt3" || $key == "dt2-dt4" ]]; then
            mask_f=~/temp_ln_subj_ttest/roi_mask/roi_dlpfc_acc_mask.nii.gz
        fi
        tcsh "$sh_basePath"/group_analysis.ttest "$proc_res_dirs" "$key" "$use_reml" "$output_dir" "$mask_f"
    }
done

# Wait for all background tasks to complete before proceeding to the next step
wait
if [[ -L ~/temp_ln_src_ttest ]]; then
    rm -rf ~/temp_ln_src_ttest # 删除创建的软链接
fi
if [[ -L ~/temp_ln_subj_ttest ]]; then
    rm -rf ~/temp_ln_subj_ttest
fi
endTime=$(date +%s)
echo "Finished in" $((endTime - beginTime)) "s" [$begin_time - $(date)]
