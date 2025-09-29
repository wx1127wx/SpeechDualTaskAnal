#!/bin/bash
# Copyright (c) 2023 Institute of Health and Medical Technology, Hefei Institutes of Physical Science, CAS, All Rights Reserved.
# @Time     : 2023/8/19
# @Author   : ZL.Z
# @Email    : zzl1124@mail.ustc.edu.cn
# @Reference: None
# @FileName : s3.batch_save_surfvol_figs.sh
# @Software : Shell; VS Code; Ubuntu 18.04.5 LTS (GNU/Linux 5.4.0-59-generic x86_64)
# @Hardware : Intel Xeon 6258R 2.7G 28C*2
# @Version  : V1.0.0 - ZL.Z: 2023/8/19
#             First version.
# @License  : None
# @Brief    : 批量保存Vol到大脑皮层的激活图

beginTime=$(date +%s)
begin_time=$(date)
# get basement path of the corrent shell file
# sh_basePath=$(cd `dirname $0`;pwd)
sh_basePath="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# 将中文路径设置软链接
if [[ -L ~/temp_ln_subj_savefig ]]; then
    rm -rf ~/temp_ln_subj_savefig
fi
ln -s "$sh_basePath" ~/temp_ln_subj_savefig # 创建软链接
# output path
output_dir=~/temp_ln_subj_savefig/results/

# anat=~/data/SUMA/suma_MNI152_2009/MNI152_2009_SurfVol.nii
anat=/home/medicalimage_admin/Applications/afni/MNI152_T1_2009c+tlrc.
spec=~/data/SUMA/suma_MNI152_2009/std.141.MNI152_2009_both.spec

# labels="group.PPI.BN.broca.bm.dlpfc_acc.nii.gz"

# running
for ppi_f in $output_dir/group.PPI.BN.*.bm.*.nii.gz; do
    {
        ppi=$(basename "$ppi_f")
        ppi_name=$(echo $ppi | sed 's/group\.PPI\.BN\.//;s/\.nii\.gz//')
        echo '----------------------------' Saving "$ppi" SurfVol Figs '----------------------------'
        fig_dir=$output_dir/SUMA_Recordings/$ppi_name
        if [[ -d "${fig_dir}" ]]; then
            rm -rf "${fig_dir}"
        fi
        mkdir -p "${fig_dir}"
        func_vol=$output_dir/$ppi
        afni -yesplugouts -niml $anat $func_vol &
        plugout_drive -com "QUIET_PLUGOUTS" \
            -com "SEE_OVERLAY +" \
            -com "SWITCH_UNDERLAY $anat" \
            -com "SWITCH_OVERLAY $ppi 24 25" \
            -com "SET_THRESHNEW 0.005 *p" \
            -quit
        sleep 1
        # 手动点击Clusterize、NN level=2、Voxels=30
        read -r -t 600 -p "等待对AFNI进行Clusterize操作,操作后请按下回车按键以继续,否则600秒后继续!"
        plugout_drive -com "QUIET_PLUGOUTS" \
            -com "PBAR_SAVEIM $output_dir/colorbar.png" \
            -com "OPEN_WINDOW axialimage mont=4x4:8" \
            -com "OPEN_WINDOW sagittalimage mont=4x4:9" \
            -com "OPEN_WINDOW coronalimage mont=4x4:9" \
            -com "SET_DICOM_XYZ -22 40 30" \
            -com "SET_XHAIRS OFF" \
            -com "SAVE_PNG axialimage $fig_dir/axial.png" \
            -com "SAVE_PNG sagittalimage $fig_dir/sagittal.png" \
            -com "SAVE_PNG coronalimage $fig_dir/coronal.png" \
            -quit
        sleep 1
        suma -spec $spec -sv $anat -niml &
        sleep 3
        # talk with afni
        DriveSuma -com viewer_cont -key t -key:r2 period -key F3 -key F4 -key F5
        sleep 3
        DriveSuma -com viewer_cont -viewer_size 1200 1000
        DriveSuma -com viewer_cont -key:d Z -key a
        sleep 5
        for pose in "Ctrl+left" "Ctrl+right" "Ctrl+up" "Ctrl+down"; do
            {
                fig_name=$(echo "$pose" | awk '{split($1, arr, "+"); print arr[2]}')
                if [[ $pose == "Ctrl+left" || $pose == "Ctrl+right" ]]; then
                    fig_f=$fig_dir/${fig_name}_lateral.jpg
                else
                    fig_f=$fig_dir/$fig_name.jpg
                fi
                if [[ $pose == "Ctrl+down" ]]; then
                    DriveSuma -com viewer_cont -key:r3:d z
                fi
                DriveSuma -com viewer_cont -key:d "$pose" -autorecord "$fig_f" -key 'Ctrl+r'
                sleep 1
            }
        done
        DriveSuma -com viewer_cont -key:d 'Ctrl+left' -key:d '[' -key:r7:d Z
        DriveSuma -com viewer_cont -autorecord $fig_dir/right_medial.jpg -key 'Ctrl+r'
        sleep 1
        DriveSuma -com viewer_cont -key:d 'Ctrl+right' -key:d '[' -key:d ']'
        DriveSuma -com viewer_cont -autorecord $fig_dir/left_medial.jpg -key 'Ctrl+r'
        sleep 3
        DriveSuma -com kill_suma
        plugout_drive -com "QUIT"
        sleep 1
        for raw_name in "$fig_dir"/*; do
            {
                if [[ ! "axial.png sagittal.png coronal.png" =~ $(basename "$raw_name") ]]; then
                    new_name=$(basename "$raw_name" | awk '{split($1, arr, ".A."); print arr[1]}')
                    mv "$raw_name" "$fig_dir/$new_name.jpg"
                fi
            }
        done
        # 获取各个任务激活对应的mask，以便在BrainNet中的标准脑上呈现
        3dcalc -a $func_vol[25] -expr "ifelse(step(a*a-2.8437*2.8437),a,0)" -prefix "$fig_dir"/"$ppi_name".nii.gz -overwrite
    }
done

# Wait for all background tasks to complete before proceeding to the next step
wait
if [[ -L ~/temp_ln_subj_savefig ]]; then
    rm -rf ~/temp_ln_subj_savefig
fi
endTime=$(date +%s)
echo "Finished in" $((endTime - beginTime)) "s" [$begin_time - $(date)]
