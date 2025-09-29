#!/bin/bash
# Copyright (c) 2023 Institute of Health and Medical Technology, Hefei Institutes of Physical Science, CAS, All Rights Reserved.
# @Time     : 2023/8/18
# @Author   : ZL.Z
# @Email    : zzl1124@mail.ustc.edu.cn
# @Reference: None
# @FileName : s1.ppi.sh
# @Software : Shell; VS Code; Ubuntu 18.04.5 LTS (GNU/Linux 5.4.0-59-generic x86_64)
# @Hardware : Intel Xeon 6258R 2.7G 28C*2
# @Version  : V1.0.0 - ZL.Z: 2023/8/18
#             First version.
# @License  : None
# @Brief    : PPI分析第一步：执行PPI分析


# note location of scripts and data
beginTime=$(date +%s)
begin_time=$(date)
# get basement path of the corrent shell file
# sh_basePath=$(cd `dirname $0`;pwd)
sh_basePath="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# 将中文路径设置软链接
# 判断主文件夹下是否存在上次设置的软链接，存在则删除
if [[ -L ~/temp_ln_proc_res ]]; then
    rm -rf ~/temp_ln_proc_res
fi
ln -s "/home/medicaldata2/ZZLData/Datasets/audio/言语与认知老化实验/言语与认知老化" ~/temp_ln_proc_res # 创建软链接
# output path
subj_dirs=~/temp_ln_proc_res/data/raw_data/t-fMRI/MRI/afni_proc_res/
# timing files path
timing_dirs=~/data/言语与认知老化实验/言语与认知老化/analysis/MRI_process/func/CreatingTimingFilesAFNI/timing_file/BLOCK/

subid_for_run="20230227087 20230308093"

# Brainnetome atlas: Broca(29,33,35,37,39);PMC(53,54,57,58,59,60,61,62,67,68);dlPFC(5,6,11,12,15,16,19,20,21,22);
# ACC(177,178,179,180,183,184,187,188);
# hippo(15,216,217,218);thalamus(231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246)
3dcalc -a roi_mask/brainnetome/BN_Atlas_246_1mm.nii.gz \
       -expr "amongst(a,29,33,35,37,39)" -prefix roi_mask/seed/BN.broca.nii.gz -overwrite
3dresample -master $subj_dirs/20220928002.results/stats.20220928002_REML+tlrc -input roi_mask/seed/BN.broca.nii.gz \
           -prefix roi_mask/seed/seed_BN.broca.nii.gz -overwrite
3dcalc -a roi_mask/brainnetome/BN_Atlas_246_1mm.nii.gz \
       -expr "amongst(a,231,233,235,237,239,241,243,245)" -prefix roi_mask/seed/BN.lh.thalamus.nii.gz -overwrite
3dresample -master $subj_dirs/20220928002.results/stats.20220928002_REML+tlrc -input roi_mask/seed/BN.lh.thalamus.nii.gz \
           -prefix roi_mask/seed/seed_BN.lh.thalamus.nii.gz -overwrite
3dcalc -a roi_mask/brainnetome/BN_Atlas_246_1mm.nii.gz \
       -expr "amongst(a,232,234,236,238,240,242,244,246)" -prefix roi_mask/seed/BN.rh.thalamus.nii.gz -overwrite
3dresample -master $subj_dirs/20220928002.results/stats.20220928002_REML+tlrc -input roi_mask/seed/BN.rh.thalamus.nii.gz \
           -prefix roi_mask/seed/seed_BN.rh.thalamus.nii.gz -overwrite

roi_seed="BN.broca BN.lh.thalamus BN.rh.thalamus"

# parallel running
for subj_dir in "$subj_dirs"/*; do
    {
        # verify existence of things
        if [[ ! -d $subj_dir ]]; then
           echo "** missing subject results directory $subj_dir"
           exit 1
        fi
        _subid=$(basename "$subj_dir")
        subid=${_subid%.results}
        # if [[ $subid_for_run =~ $subid ]]; then
        if [[ ! -d $subj_dir/stimuli/BLOCK ]]; then
            mkdir $subj_dir/stimuli/BLOCK
            cp $timing_dirs/"$subid"/timing_*.txt $subj_dir/stimuli/BLOCK
        fi
        # do all of the work in the $subid.results directory...
        cd $subj_dir
        for seed in $roi_seed; do  # start with seed at ROI
            {
                echo '----------------------------' "$subid": $seed '----------------------------'
                # echo 4 28 16 | 3dUndump -xyz -srad 5 -master stats."$subid"_REML+tlrc -prefix mask.$seed.nii.gz -
                # generate seed time series, ppi.seed.1D (note that mask dset is unneeded, but visually useful)
                3dmaskave -quiet -mask $sh_basePath/roi_mask/seed/seed_"$seed".nii.gz errts."$subid"_REML+tlrc > ppi.seed.$seed.1D
                # generate PPI regressors from seed and timing files
                tcsh $sh_basePath/cmd.ppi.2.make.regs $seed
                # and copy the results into the stimuli directory
                cp work.$seed/p6.* ppi.seed.$seed.1D stimuli
                # and just to see consider:
                #    1dplot -one ppi.seed.$seed.1D work.$seed/p7.$seed.sum.PPI.1D
                #    1dplot ppi.seed.$seed.1D work.$seed/p6.*
                # run a 3dDeconvolve command for the PPI
                tcsh $sh_basePath/proc.3dd.ppi.post.full $subid $seed
            } &
        done
        wait
        # fi
    }
done

# Wait for all background tasks to complete before proceeding to the next step
wait
if [[ -L ~/temp_ln_proc_res ]]; then
    rm -rf ~/temp_ln_proc_res # 删除创建的软链接
fi
endTime=$(date +%s)
echo "Finished in" $((endTime - beginTime)) "s" [$begin_time - $(date)]



