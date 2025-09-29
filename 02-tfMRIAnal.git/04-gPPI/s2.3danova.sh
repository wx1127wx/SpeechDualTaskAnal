#!/bin/bash
# Copyright (c) 2023 Institute of Health and Medical Technology, Hefei Institutes of Physical Science, CAS, All Rights Reserved.
# @Time     : 2023/8/19
# @Author   : ZL.Z
# @Email    : zzl1124@mail.ustc.edu.cn
# @Reference: None
# @FileName : s2.3danova2.sh
# @Software : Shell; VS Code; Ubuntu 18.04.5 LTS (GNU/Linux 5.4.0-59-generic x86_64)
# @Hardware : Intel Xeon 6258R 2.7G 28C*2
# @Version  : V1.0.0 - ZL.Z: 2023/8/19
#             First version.
# @License  : None
# @Brief    : PPI组分析：ANOVA

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

if [[ -L ~/temp_ln_subj_anova ]]; then
    rm -rf ~/temp_ln_subj_anova
fi
ln -s "$sh_basePath" ~/temp_ln_subj_anova # 创建软链接
# output path
output_dir=~/temp_ln_subj_anova/results/
find $output_dir -type f -name "group.PPI.*.nii.gz" -exec rm {} \;
# Brainnetome atlas: Broca(29,33,35,37,39);dlPFC(5,6,11,12,15,16,19,20,21,22);ACC(177,178,179,180,183,184,187,188);
# thalamus(231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246)
3dcalc -a roi_mask/brainnetome/BN_Atlas_246_1mm.nii.gz \
       -expr "amongst(a,29,33,35,37,39)" -prefix roi_mask/roi_BN.broca.nii.gz -overwrite
3dresample -master $subj_dirs/20220928002.results/stats.20220928002_REML+tlrc -input roi_mask/roi_BN.broca.nii.gz \
           -prefix roi_mask/roi_mask_BN.broca.nii.gz -overwrite
3dcalc -a roi_mask/brainnetome/BN_Atlas_246_1mm.nii.gz \
       -expr "amongst(a,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246)" \
       -prefix roi_mask/roi_BN.thalamus.nii.gz -overwrite
3dresample -master $subj_dirs/20220928002.results/stats.20220928002_REML+tlrc -input roi_mask/roi_BN.thalamus.nii.gz \
           -prefix roi_mask/roi_mask_BN.thalamus.nii.gz -overwrite
3dcalc -a roi_mask/brainnetome/BN_Atlas_246_1mm.nii.gz \
       -expr "amongst(a,5,6,11,12,15,16,19,20,21,22,177,178,179,180,183,184,187,188)" \
       -prefix roi_mask/roi_BN.dlpfc_acc.nii.gz -overwrite
3dresample -master $subj_dirs/20220928002.results/stats.20220928002_REML+tlrc -input roi_mask/roi_BN.dlpfc_acc.nii.gz \
           -prefix roi_mask/roi_mask_BN.dlpfc_acc.nii.gz -overwrite

roi_seed="BN.broca BN.lh.thalamus BN.rh.thalamus"

for seed in $roi_seed; do
    {
        echo '----------------------------' ANOVA: "$seed" '----------------------------'
        # one factor: task type (5 levels)
        3dANOVA -levels 5                                                                            \
-dset 1 "$subj_dirs/20220928003.results/PPI.$seed.stats.20220928003+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20220928003.results/PPI.$seed.stats.20220928003+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20220928003.results/PPI.$seed.stats.20220928003+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20220928003.results/PPI.$seed.stats.20220928003+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20220928003.results/PPI.$seed.stats.20220928003+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20220929007.results/PPI.$seed.stats.20220929007+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20220929007.results/PPI.$seed.stats.20220929007+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20220929007.results/PPI.$seed.stats.20220929007+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20220929007.results/PPI.$seed.stats.20220929007+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20220929007.results/PPI.$seed.stats.20220929007+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20220928002.results/PPI.$seed.stats.20220928002+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20220928002.results/PPI.$seed.stats.20220928002+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20220928002.results/PPI.$seed.stats.20220928002+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20220928002.results/PPI.$seed.stats.20220928002+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20220928002.results/PPI.$seed.stats.20220928002+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20221106032.results/PPI.$seed.stats.20221106032+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20221106032.results/PPI.$seed.stats.20221106032+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20221106032.results/PPI.$seed.stats.20221106032+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20221106032.results/PPI.$seed.stats.20221106032+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20221106032.results/PPI.$seed.stats.20221106032+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20221023023.results/PPI.$seed.stats.20221023023+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20221023023.results/PPI.$seed.stats.20221023023+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20221023023.results/PPI.$seed.stats.20221023023+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20221023023.results/PPI.$seed.stats.20221023023+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20221023023.results/PPI.$seed.stats.20221023023+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20221028026.results/PPI.$seed.stats.20221028026+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20221028026.results/PPI.$seed.stats.20221028026+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20221028026.results/PPI.$seed.stats.20221028026+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20221028026.results/PPI.$seed.stats.20221028026+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20221028026.results/PPI.$seed.stats.20221028026+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20221105031.results/PPI.$seed.stats.20221105031+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20221105031.results/PPI.$seed.stats.20221105031+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20221105031.results/PPI.$seed.stats.20221105031+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20221105031.results/PPI.$seed.stats.20221105031+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20221105031.results/PPI.$seed.stats.20221105031+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20221005018.results/PPI.$seed.stats.20221005018+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20221005018.results/PPI.$seed.stats.20221005018+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20221005018.results/PPI.$seed.stats.20221005018+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20221005018.results/PPI.$seed.stats.20221005018+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20221005018.results/PPI.$seed.stats.20221005018+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20221205036.results/PPI.$seed.stats.20221205036+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20221205036.results/PPI.$seed.stats.20221205036+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20221205036.results/PPI.$seed.stats.20221205036+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20221205036.results/PPI.$seed.stats.20221205036+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20221205036.results/PPI.$seed.stats.20221205036+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20221022022.results/PPI.$seed.stats.20221022022+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20221022022.results/PPI.$seed.stats.20221022022+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20221022022.results/PPI.$seed.stats.20221022022+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20221022022.results/PPI.$seed.stats.20221022022+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20221022022.results/PPI.$seed.stats.20221022022+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20221206039.results/PPI.$seed.stats.20221206039+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20221206039.results/PPI.$seed.stats.20221206039+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20221206039.results/PPI.$seed.stats.20221206039+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20221206039.results/PPI.$seed.stats.20221206039+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20221206039.results/PPI.$seed.stats.20221206039+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20221206038.results/PPI.$seed.stats.20221206038+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20221206038.results/PPI.$seed.stats.20221206038+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20221206038.results/PPI.$seed.stats.20221206038+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20221206038.results/PPI.$seed.stats.20221206038+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20221206038.results/PPI.$seed.stats.20221206038+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20221001015.results/PPI.$seed.stats.20221001015+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20221001015.results/PPI.$seed.stats.20221001015+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20221001015.results/PPI.$seed.stats.20221001015+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20221001015.results/PPI.$seed.stats.20221001015+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20221001015.results/PPI.$seed.stats.20221001015+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20221104029.results/PPI.$seed.stats.20221104029+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20221104029.results/PPI.$seed.stats.20221104029+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20221104029.results/PPI.$seed.stats.20221104029+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20221104029.results/PPI.$seed.stats.20221104029+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20221104029.results/PPI.$seed.stats.20221104029+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20221023024.results/PPI.$seed.stats.20221023024+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20221023024.results/PPI.$seed.stats.20221023024+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20221023024.results/PPI.$seed.stats.20221023024+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20221023024.results/PPI.$seed.stats.20221023024+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20221023024.results/PPI.$seed.stats.20221023024+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20221204034.results/PPI.$seed.stats.20221204034+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20221204034.results/PPI.$seed.stats.20221204034+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20221204034.results/PPI.$seed.stats.20221204034+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20221204034.results/PPI.$seed.stats.20221204034+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20221204034.results/PPI.$seed.stats.20221204034+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230208049.results/PPI.$seed.stats.20230208049+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230208049.results/PPI.$seed.stats.20230208049+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230208049.results/PPI.$seed.stats.20230208049+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230208049.results/PPI.$seed.stats.20230208049+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230208049.results/PPI.$seed.stats.20230208049+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20220929006.results/PPI.$seed.stats.20220929006+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20220929006.results/PPI.$seed.stats.20220929006+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20220929006.results/PPI.$seed.stats.20220929006+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20220929006.results/PPI.$seed.stats.20220929006+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20220929006.results/PPI.$seed.stats.20220929006+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230226084.results/PPI.$seed.stats.20230226084+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230226084.results/PPI.$seed.stats.20230226084+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230226084.results/PPI.$seed.stats.20230226084+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230226084.results/PPI.$seed.stats.20230226084+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230226084.results/PPI.$seed.stats.20230226084+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20221001016.results/PPI.$seed.stats.20221001016+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20221001016.results/PPI.$seed.stats.20221001016+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20221001016.results/PPI.$seed.stats.20221001016+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20221001016.results/PPI.$seed.stats.20221001016+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20221001016.results/PPI.$seed.stats.20221001016+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230214070.results/PPI.$seed.stats.20230214070+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230214070.results/PPI.$seed.stats.20230214070+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230214070.results/PPI.$seed.stats.20230214070+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230214070.results/PPI.$seed.stats.20230214070+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230214070.results/PPI.$seed.stats.20230214070+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20221025025.results/PPI.$seed.stats.20221025025+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20221025025.results/PPI.$seed.stats.20221025025+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20221025025.results/PPI.$seed.stats.20221025025+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20221025025.results/PPI.$seed.stats.20221025025+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20221025025.results/PPI.$seed.stats.20221025025+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230211056.results/PPI.$seed.stats.20230211056+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230211056.results/PPI.$seed.stats.20230211056+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230211056.results/PPI.$seed.stats.20230211056+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230211056.results/PPI.$seed.stats.20230211056+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230211056.results/PPI.$seed.stats.20230211056+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20221020021.results/PPI.$seed.stats.20221020021+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20221020021.results/PPI.$seed.stats.20221020021+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20221020021.results/PPI.$seed.stats.20221020021+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20221020021.results/PPI.$seed.stats.20221020021+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20221020021.results/PPI.$seed.stats.20221020021+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230211058.results/PPI.$seed.stats.20230211058+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230211058.results/PPI.$seed.stats.20230211058+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230211058.results/PPI.$seed.stats.20230211058+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230211058.results/PPI.$seed.stats.20230211058+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230211058.results/PPI.$seed.stats.20230211058+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20221017020.results/PPI.$seed.stats.20221017020+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20221017020.results/PPI.$seed.stats.20221017020+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20221017020.results/PPI.$seed.stats.20221017020+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20221017020.results/PPI.$seed.stats.20221017020+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20221017020.results/PPI.$seed.stats.20221017020+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230211059.results/PPI.$seed.stats.20230211059+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230211059.results/PPI.$seed.stats.20230211059+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230211059.results/PPI.$seed.stats.20230211059+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230211059.results/PPI.$seed.stats.20230211059+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230211059.results/PPI.$seed.stats.20230211059+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20221110033.results/PPI.$seed.stats.20221110033+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20221110033.results/PPI.$seed.stats.20221110033+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20221110033.results/PPI.$seed.stats.20221110033+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20221110033.results/PPI.$seed.stats.20221110033+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20221110033.results/PPI.$seed.stats.20221110033+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20221016019.results/PPI.$seed.stats.20221016019+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20221016019.results/PPI.$seed.stats.20221016019+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20221016019.results/PPI.$seed.stats.20221016019+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20221016019.results/PPI.$seed.stats.20221016019+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20221016019.results/PPI.$seed.stats.20221016019+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230211057.results/PPI.$seed.stats.20230211057+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230211057.results/PPI.$seed.stats.20230211057+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230211057.results/PPI.$seed.stats.20230211057+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230211057.results/PPI.$seed.stats.20230211057+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230211057.results/PPI.$seed.stats.20230211057+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230214069.results/PPI.$seed.stats.20230214069+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230214069.results/PPI.$seed.stats.20230214069+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230214069.results/PPI.$seed.stats.20230214069+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230214069.results/PPI.$seed.stats.20230214069+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230214069.results/PPI.$seed.stats.20230214069+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20221205037.results/PPI.$seed.stats.20221205037+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20221205037.results/PPI.$seed.stats.20221205037+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20221205037.results/PPI.$seed.stats.20221205037+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20221205037.results/PPI.$seed.stats.20221205037+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20221205037.results/PPI.$seed.stats.20221205037+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230108045.results/PPI.$seed.stats.20230108045+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230108045.results/PPI.$seed.stats.20230108045+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230108045.results/PPI.$seed.stats.20230108045+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230108045.results/PPI.$seed.stats.20230108045+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230108045.results/PPI.$seed.stats.20230108045+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20221102030.results/PPI.$seed.stats.20221102030+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20221102030.results/PPI.$seed.stats.20221102030+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20221102030.results/PPI.$seed.stats.20221102030+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20221102030.results/PPI.$seed.stats.20221102030+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20221102030.results/PPI.$seed.stats.20221102030+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20221213042.results/PPI.$seed.stats.20221213042+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20221213042.results/PPI.$seed.stats.20221213042+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20221213042.results/PPI.$seed.stats.20221213042+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20221213042.results/PPI.$seed.stats.20221213042+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20221213042.results/PPI.$seed.stats.20221213042+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230208051.results/PPI.$seed.stats.20230208051+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230208051.results/PPI.$seed.stats.20230208051+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230208051.results/PPI.$seed.stats.20230208051+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230208051.results/PPI.$seed.stats.20230208051+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230208051.results/PPI.$seed.stats.20230208051+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230208050.results/PPI.$seed.stats.20230208050+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230208050.results/PPI.$seed.stats.20230208050+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230208050.results/PPI.$seed.stats.20230208050+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230208050.results/PPI.$seed.stats.20230208050+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230208050.results/PPI.$seed.stats.20230208050+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230216071.results/PPI.$seed.stats.20230216071+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230216071.results/PPI.$seed.stats.20230216071+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230216071.results/PPI.$seed.stats.20230216071+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230216071.results/PPI.$seed.stats.20230216071+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230216071.results/PPI.$seed.stats.20230216071+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20220929010.results/PPI.$seed.stats.20220929010+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20220929010.results/PPI.$seed.stats.20220929010+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20220929010.results/PPI.$seed.stats.20220929010+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20220929010.results/PPI.$seed.stats.20220929010+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20220929010.results/PPI.$seed.stats.20220929010+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230224082.results/PPI.$seed.stats.20230224082+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230224082.results/PPI.$seed.stats.20230224082+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230224082.results/PPI.$seed.stats.20230224082+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230224082.results/PPI.$seed.stats.20230224082+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230224082.results/PPI.$seed.stats.20230224082+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20221205035.results/PPI.$seed.stats.20221205035+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20221205035.results/PPI.$seed.stats.20221205035+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20221205035.results/PPI.$seed.stats.20221205035+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20221205035.results/PPI.$seed.stats.20221205035+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20221205035.results/PPI.$seed.stats.20221205035+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20221028027.results/PPI.$seed.stats.20221028027+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20221028027.results/PPI.$seed.stats.20221028027+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20221028027.results/PPI.$seed.stats.20221028027+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20221028027.results/PPI.$seed.stats.20221028027+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20221028027.results/PPI.$seed.stats.20221028027+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230307092.results/PPI.$seed.stats.20230307092+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230307092.results/PPI.$seed.stats.20230307092+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230307092.results/PPI.$seed.stats.20230307092+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230307092.results/PPI.$seed.stats.20230307092+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230307092.results/PPI.$seed.stats.20230307092+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230220075.results/PPI.$seed.stats.20230220075+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230220075.results/PPI.$seed.stats.20230220075+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230220075.results/PPI.$seed.stats.20230220075+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230220075.results/PPI.$seed.stats.20230220075+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230220075.results/PPI.$seed.stats.20230220075+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230212061.results/PPI.$seed.stats.20230212061+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230212061.results/PPI.$seed.stats.20230212061+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230212061.results/PPI.$seed.stats.20230212061+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230212061.results/PPI.$seed.stats.20230212061+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230212061.results/PPI.$seed.stats.20230212061+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230213064.results/PPI.$seed.stats.20230213064+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230213064.results/PPI.$seed.stats.20230213064+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230213064.results/PPI.$seed.stats.20230213064+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230213064.results/PPI.$seed.stats.20230213064+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230213064.results/PPI.$seed.stats.20230213064+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230213066.results/PPI.$seed.stats.20230213066+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230213066.results/PPI.$seed.stats.20230213066+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230213066.results/PPI.$seed.stats.20230213066+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230213066.results/PPI.$seed.stats.20230213066+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230213066.results/PPI.$seed.stats.20230213066+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230222079.results/PPI.$seed.stats.20230222079+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230222079.results/PPI.$seed.stats.20230222079+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230222079.results/PPI.$seed.stats.20230222079+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230222079.results/PPI.$seed.stats.20230222079+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230222079.results/PPI.$seed.stats.20230222079+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230213067.results/PPI.$seed.stats.20230213067+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230213067.results/PPI.$seed.stats.20230213067+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230213067.results/PPI.$seed.stats.20230213067+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230213067.results/PPI.$seed.stats.20230213067+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230213067.results/PPI.$seed.stats.20230213067+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230108044.results/PPI.$seed.stats.20230108044+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230108044.results/PPI.$seed.stats.20230108044+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230108044.results/PPI.$seed.stats.20230108044+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230108044.results/PPI.$seed.stats.20230108044+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230108044.results/PPI.$seed.stats.20230108044+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230315098.results/PPI.$seed.stats.20230315098+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230315098.results/PPI.$seed.stats.20230315098+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230315098.results/PPI.$seed.stats.20230315098+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230315098.results/PPI.$seed.stats.20230315098+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230315098.results/PPI.$seed.stats.20230315098+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230211055.results/PPI.$seed.stats.20230211055+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230211055.results/PPI.$seed.stats.20230211055+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230211055.results/PPI.$seed.stats.20230211055+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230211055.results/PPI.$seed.stats.20230211055+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230211055.results/PPI.$seed.stats.20230211055+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230211060.results/PPI.$seed.stats.20230211060+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230211060.results/PPI.$seed.stats.20230211060+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230211060.results/PPI.$seed.stats.20230211060+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230211060.results/PPI.$seed.stats.20230211060+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230211060.results/PPI.$seed.stats.20230211060+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230303091.results/PPI.$seed.stats.20230303091+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230303091.results/PPI.$seed.stats.20230303091+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230303091.results/PPI.$seed.stats.20230303091+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230303091.results/PPI.$seed.stats.20230303091+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230303091.results/PPI.$seed.stats.20230303091+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230209054.results/PPI.$seed.stats.20230209054+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230209054.results/PPI.$seed.stats.20230209054+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230209054.results/PPI.$seed.stats.20230209054+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230209054.results/PPI.$seed.stats.20230209054+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230209054.results/PPI.$seed.stats.20230209054+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20220929009.results/PPI.$seed.stats.20220929009+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20220929009.results/PPI.$seed.stats.20220929009+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20220929009.results/PPI.$seed.stats.20220929009+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20220929009.results/PPI.$seed.stats.20220929009+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20220929009.results/PPI.$seed.stats.20220929009+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20221029028.results/PPI.$seed.stats.20221029028+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20221029028.results/PPI.$seed.stats.20221029028+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20221029028.results/PPI.$seed.stats.20221029028+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20221029028.results/PPI.$seed.stats.20221029028+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20221029028.results/PPI.$seed.stats.20221029028+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230108046.results/PPI.$seed.stats.20230108046+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230108046.results/PPI.$seed.stats.20230108046+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230108046.results/PPI.$seed.stats.20230108046+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230108046.results/PPI.$seed.stats.20230108046+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230108046.results/PPI.$seed.stats.20230108046+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230213065.results/PPI.$seed.stats.20230213065+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230213065.results/PPI.$seed.stats.20230213065+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230213065.results/PPI.$seed.stats.20230213065+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230213065.results/PPI.$seed.stats.20230213065+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230213065.results/PPI.$seed.stats.20230213065+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230212062.results/PPI.$seed.stats.20230212062+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230212062.results/PPI.$seed.stats.20230212062+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230212062.results/PPI.$seed.stats.20230212062+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230212062.results/PPI.$seed.stats.20230212062+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230212062.results/PPI.$seed.stats.20230212062+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230209053.results/PPI.$seed.stats.20230209053+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230209053.results/PPI.$seed.stats.20230209053+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230209053.results/PPI.$seed.stats.20230209053+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230209053.results/PPI.$seed.stats.20230209053+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230209053.results/PPI.$seed.stats.20230209053+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230312094.results/PPI.$seed.stats.20230312094+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230312094.results/PPI.$seed.stats.20230312094+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230312094.results/PPI.$seed.stats.20230312094+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230312094.results/PPI.$seed.stats.20230312094+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230312094.results/PPI.$seed.stats.20230312094+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20221213041.results/PPI.$seed.stats.20221213041+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20221213041.results/PPI.$seed.stats.20221213041+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20221213041.results/PPI.$seed.stats.20221213041+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20221213041.results/PPI.$seed.stats.20221213041+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20221213041.results/PPI.$seed.stats.20221213041+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230208047.results/PPI.$seed.stats.20230208047+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230208047.results/PPI.$seed.stats.20230208047+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230208047.results/PPI.$seed.stats.20230208047+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230208047.results/PPI.$seed.stats.20230208047+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230208047.results/PPI.$seed.stats.20230208047+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230219073.results/PPI.$seed.stats.20230219073+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230219073.results/PPI.$seed.stats.20230219073+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230219073.results/PPI.$seed.stats.20230219073+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230219073.results/PPI.$seed.stats.20230219073+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230219073.results/PPI.$seed.stats.20230219073+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230316100.results/PPI.$seed.stats.20230316100+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230316100.results/PPI.$seed.stats.20230316100+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230316100.results/PPI.$seed.stats.20230316100+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230316100.results/PPI.$seed.stats.20230316100+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230316100.results/PPI.$seed.stats.20230316100+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230208048.results/PPI.$seed.stats.20230208048+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230208048.results/PPI.$seed.stats.20230208048+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230208048.results/PPI.$seed.stats.20230208048+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230208048.results/PPI.$seed.stats.20230208048+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230208048.results/PPI.$seed.stats.20230208048+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230216072.results/PPI.$seed.stats.20230216072+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230216072.results/PPI.$seed.stats.20230216072+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230216072.results/PPI.$seed.stats.20230216072+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230216072.results/PPI.$seed.stats.20230216072+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230216072.results/PPI.$seed.stats.20230216072+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20221019017.results/PPI.$seed.stats.20221019017+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20221019017.results/PPI.$seed.stats.20221019017+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20221019017.results/PPI.$seed.stats.20221019017+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20221019017.results/PPI.$seed.stats.20221019017+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20221019017.results/PPI.$seed.stats.20221019017+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230227088.results/PPI.$seed.stats.20230227088+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230227088.results/PPI.$seed.stats.20230227088+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230227088.results/PPI.$seed.stats.20230227088+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230227088.results/PPI.$seed.stats.20230227088+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230227088.results/PPI.$seed.stats.20230227088+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230226083.results/PPI.$seed.stats.20230226083+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230226083.results/PPI.$seed.stats.20230226083+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230226083.results/PPI.$seed.stats.20230226083+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230226083.results/PPI.$seed.stats.20230226083+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230226083.results/PPI.$seed.stats.20230226083+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230227087.results/PPI.$seed.stats.20230227087+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230227087.results/PPI.$seed.stats.20230227087+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230227087.results/PPI.$seed.stats.20230227087+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230227087.results/PPI.$seed.stats.20230227087+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230227087.results/PPI.$seed.stats.20230227087+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20221209040.results/PPI.$seed.stats.20221209040+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20221209040.results/PPI.$seed.stats.20221209040+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20221209040.results/PPI.$seed.stats.20221209040+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20221209040.results/PPI.$seed.stats.20221209040+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20221209040.results/PPI.$seed.stats.20221209040+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230314097.results/PPI.$seed.stats.20230314097+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230314097.results/PPI.$seed.stats.20230314097+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230314097.results/PPI.$seed.stats.20230314097+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230314097.results/PPI.$seed.stats.20230314097+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230314097.results/PPI.$seed.stats.20230314097+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230223081.results/PPI.$seed.stats.20230223081+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230223081.results/PPI.$seed.stats.20230223081+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230223081.results/PPI.$seed.stats.20230223081+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230223081.results/PPI.$seed.stats.20230223081+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230223081.results/PPI.$seed.stats.20230223081+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230221078.results/PPI.$seed.stats.20230221078+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230221078.results/PPI.$seed.stats.20230221078+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230221078.results/PPI.$seed.stats.20230221078+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230221078.results/PPI.$seed.stats.20230221078+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230221078.results/PPI.$seed.stats.20230221078+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230312095.results/PPI.$seed.stats.20230312095+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230312095.results/PPI.$seed.stats.20230312095+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230312095.results/PPI.$seed.stats.20230312095+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230312095.results/PPI.$seed.stats.20230312095+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230312095.results/PPI.$seed.stats.20230312095+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230424101.results/PPI.$seed.stats.20230424101+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230424101.results/PPI.$seed.stats.20230424101+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230424101.results/PPI.$seed.stats.20230424101+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230424101.results/PPI.$seed.stats.20230424101+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230424101.results/PPI.$seed.stats.20230424101+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230214068.results/PPI.$seed.stats.20230214068+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230214068.results/PPI.$seed.stats.20230214068+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230214068.results/PPI.$seed.stats.20230214068+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230214068.results/PPI.$seed.stats.20230214068+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230214068.results/PPI.$seed.stats.20230214068+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230212063.results/PPI.$seed.stats.20230212063+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230212063.results/PPI.$seed.stats.20230212063+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230212063.results/PPI.$seed.stats.20230212063+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230212063.results/PPI.$seed.stats.20230212063+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230212063.results/PPI.$seed.stats.20230212063+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230424102.results/PPI.$seed.stats.20230424102+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230424102.results/PPI.$seed.stats.20230424102+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230424102.results/PPI.$seed.stats.20230424102+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230424102.results/PPI.$seed.stats.20230424102+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230424102.results/PPI.$seed.stats.20230424102+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230315099.results/PPI.$seed.stats.20230315099+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230315099.results/PPI.$seed.stats.20230315099+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230315099.results/PPI.$seed.stats.20230315099+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230315099.results/PPI.$seed.stats.20230315099+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230315099.results/PPI.$seed.stats.20230315099+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20221214043.results/PPI.$seed.stats.20221214043+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20221214043.results/PPI.$seed.stats.20221214043+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20221214043.results/PPI.$seed.stats.20221214043+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20221214043.results/PPI.$seed.stats.20221214043+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20221214043.results/PPI.$seed.stats.20221214043+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230228089.results/PPI.$seed.stats.20230228089+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230228089.results/PPI.$seed.stats.20230228089+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230228089.results/PPI.$seed.stats.20230228089+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230228089.results/PPI.$seed.stats.20230228089+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230228089.results/PPI.$seed.stats.20230228089+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230308093.results/PPI.$seed.stats.20230308093+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230308093.results/PPI.$seed.stats.20230308093+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230308093.results/PPI.$seed.stats.20230308093+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230308093.results/PPI.$seed.stats.20230308093+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230308093.results/PPI.$seed.stats.20230308093+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230220076.results/PPI.$seed.stats.20230220076+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230220076.results/PPI.$seed.stats.20230220076+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230220076.results/PPI.$seed.stats.20230220076+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230220076.results/PPI.$seed.stats.20230220076+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230220076.results/PPI.$seed.stats.20230220076+tlrc.HEAD[PPI.dt4#0_Coef]"   \
-dset 1 "$subj_dirs/20230226086.results/PPI.$seed.stats.20230226086+tlrc.HEAD[PPI.base#0_Coef]"   \
-dset 2 "$subj_dirs/20230226086.results/PPI.$seed.stats.20230226086+tlrc.HEAD[PPI.dt1#0_Coef]"   \
-dset 3 "$subj_dirs/20230226086.results/PPI.$seed.stats.20230226086+tlrc.HEAD[PPI.dt2#0_Coef]"   \
-dset 4 "$subj_dirs/20230226086.results/PPI.$seed.stats.20230226086+tlrc.HEAD[PPI.dt3#0_Coef]"   \
-dset 5 "$subj_dirs/20230226086.results/PPI.$seed.stats.20230226086+tlrc.HEAD[PPI.dt4#0_Coef]"   \
            -ftr Tasks -mean 1 PPI_mean:B -mean 2 PPI_mean:D1 -mean 3 PPI_mean:D2                \
            -mean 4 PPI_mean:D3 -mean 5 PPI_mean:D4 -diff 2 1 PPI_diff:D1vsB                     \
            -diff 3 1 PPI_diff:D2vsB -diff 2 4 PPI_diff:D1vsD3 -diff 3 5 PPI_diff:D2vsD4         \
            -contr -0.5 1 0 -0.5 0 PPI_contr:D1vsBD3 -contr -0.5 0 1 0 -0.5 PPI_contr:D2vsBD4    \
            -contr 0 0.5 0.5 0 0 PPI_contr:D1D2vs0 -contr -1 0.5 0.5 0 0 PPI_contr:D1D2vsB       \
            -contr 0 0.5 0.5 -0.5 -0.5 PPI_contr:D1D2vsD3D4 -contr -1 1.5 1.5 -1 -1 PPI_contr:D1D2vsBD3D4      \
            -bucket $output_dir/group.PPI.$seed.nii.gz
    sd=${seed##*.}
    for roimask in roi_mask/roi_mask*.nii.gz; do
        {
            roi=$(echo $roimask | cut -d'.' -f2)
            if [[ $sd != $roi ]]; then
               3dcalc -a $output_dir/group.PPI.$seed.nii.gz -b $roimask \
                      -expr "a*step(b)" -prefix $output_dir/group.PPI.$seed.bm.$roi.nii.gz -overwrite
            fi
        }
       done
    }
done

# Wait for all background tasks to complete before proceeding to the next step
wait
if [[ -L ~/temp_ln_proc_res ]]; then
    rm -rf ~/temp_ln_proc_res # 删除创建的软链接
fi
if [[ -L ~/temp_ln_subj_anova ]]; then
    rm -rf ~/temp_ln_subj_anova
fi
endTime=$(date +%s)
echo "Finished in" $((endTime - beginTime)) "s" [$begin_time - $(date)]
