#!/bin/bash
# Copyright (c) 2023 Institute of Health and Medical Technology, Hefei Institutes of Physical Science, CAS, All Rights Reserved.
# @Time     : 2023/6/30 15:00
# @Author   : ZL.Z
# @Email    : zzl1124@mail.ustc.edu.cn
# @Reference: None
# @FileName : afni_proc.sh
# @Software : Shell; VS Code; Ubuntu 18.04.5 LTS (GNU/Linux 5.4.0-59-generic x86_64)
# @Hardware : Intel Xeon 6258R 2.7G 28C*2
# @Version  : V1.0.0 - ZL.Z: 2023/6/30
#             First version.
# @License  : None
# @Brief    : 利用afni_proc.py生成AFNI的预处理脚本文件(以subj_exp示例被试生成proc.mvpa，需要修改目录才能进行后续的批处理)

beginTime=$(date +%s)
# get basement path of the corrent shell file
# sh_basePath=$(cd `dirname $0`;pwd)
sh_basePath="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

subj_exp="$sh_basePath"/example

afni_proc.py -subj_id $subj_exp -script proc.default_subid                         \
    -scr_overwrite -blocks tshift align tlrc volreg blur mask scale regress            \
    -copy_anat ${subj_exp}/anat/T1W_3D_MPRAGE_1mm_SENSE.nii.gz \
    -dsets ${subj_exp}/func/DDKDT_run*_SENSE.nii.gz            \
    -tcat_remove_first_trs 0 -align_opts_aea -giant_move       \
    -tlrc_base MNI_avg152T1+tlrc -tlrc_opts_at -init_xform AUTO_CENTER -volreg_align_to MIN_OUTLIER  \
    -volreg_align_e2a -volreg_tlrc_warp -blur_size 4.0         \
    -regress_stim_times ${subj_exp}/timing/timing_*.txt        \
    -regress_stim_labels base dt1 dt2 dt3 dt4                  \
    -regress_basis dmBLOCK -regress_stim_types AM2 -regress_censor_motion 2.0          \
    -regress_motion_per_run -regress_opts_3dD -jobs 4 -gltsym 'SYM: dt1 -base'         \
    -glt_label 1 dt1-base -gltsym 'SYM: dt2 -base' -glt_label 2 dt2-base               \
    -gltsym 'SYM: dt1 -dt3' -glt_label 3 dt1-dt3 -gltsym 'SYM: dt2 -dt4' -glt_label 4 dt2-dt4 \
    -gltsym 'SYM: dt1 -base -dt3' -glt_label 5 dt1-base-dt3 -gltsym 'SYM: dt2 -base -dt4' -glt_label 6 dt2-base-dt4 \
    -gltsym 'SYM: dt1 -dt3 -base' -glt_label 7 dt1-dt3-base -gltsym 'SYM: dt2 -dt4 -base' -glt_label 8 dt2-dt4-base \
    -regress_reml_exec -regress_make_ideal_sum sum_ideal.1D                            \
    -regress_est_blur_epits -regress_est_blur_errts

# afni_proc.py -subj_id $subj_exp -script proc.default_subid                         \
#     -scr_overwrite -blocks tshift align tlrc volreg blur mask scale regress            \
#     -copy_anat ${subj_exp}/anat/T1W_3D_MPRAGE_1mm_SENSE.nii.gz \
#     -dsets ${subj_exp}/func/DDKDT_run*_SENSE.nii.gz            \
#     -tcat_remove_first_trs 0 -align_opts_aea -giant_move       \
#     -tlrc_base MNI_avg152T1+tlrc -tlrc_opts_at -init_xform AUTO_CENTER -volreg_align_to MIN_OUTLIER  \
#     -volreg_align_e2a -volreg_tlrc_warp -blur_size 4.0         \
#     -regress_stim_times ${subj_exp}/timing/timing_*.txt        \
#     -regress_stim_labels base dt1 dt2 dt3 dt4                  \
#     -regress_basis 'BLOCK(3.2,1)' -regress_censor_motion 2.0          \
#     -regress_motion_per_run -regress_opts_3dD -jobs 4 -gltsym 'SYM: dt1 -base'         \
#     -glt_label 1 dt1-base -gltsym 'SYM: dt2 -base' -glt_label 2 dt2-base               \
#     -gltsym 'SYM: dt1 -dt3' -glt_label 3 dt1-dt3 -gltsym 'SYM: dt2 -dt4' -glt_label 4 dt2-dt4 \
#     -gltsym 'SYM: dt1 -base -dt3' -glt_label 5 dt1-base-dt3 -gltsym 'SYM: dt2 -base -dt4' -glt_label 6 dt2-base-dt4 \
#     -gltsym 'SYM: dt1 -dt3 -base' -glt_label 7 dt1-dt3-base -gltsym 'SYM: dt2 -dt4 -base' -glt_label 8 dt2-dt4-base \
#     -regress_reml_exec -regress_make_ideal_sum sum_ideal.1D                            \
#     -regress_est_blur_epits -regress_est_blur_errts #-volreg_warp_dxyz 2.0

endTime=$(date +%s)
echo "Finished in" $((endTime - beginTime)) "s"
