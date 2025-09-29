#!/bin/bash
# Copyright (c) 2022 Institute of Health and Medical Technology, Hefei Institutes of Physical Science, CAS, All Rights Reserved.
# @Time     : 2022/12/23 18:17:13
# @Author   : ZL.Z
# @Email    : zzl1124@mail.ustc.edu.cn
# @Reference: None
# @FileName : batch_save_surfvol_figs.sh
# @Software : Shell; VS Code; Ubuntu 18.04.5 LTS (GNU/Linux 5.4.0-59-generic x86_64)
# @Hardware : Intel Xeon 6258R 2.7G 28C*2
# @Version  : V1.1.0 - ZL.Z: 2023/5/25
#             增加保存Vol激活图
#             V1.0.0 - ZL.Z: 2022/12/23
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

anat=~/data/SUMA/suma_MNI152_2009/MNI152_2009_SurfVol.nii
spec=~/data/SUMA/suma_MNI152_2009/std.141.MNI152_2009_both.spec

# labels="base dt1 dt2 dt3 dt4 dt1-base dt2-base dt1-dt3 dt2-dt4 dt1-base-dt3 dt2-base-dt4 dt1-dt3-base dt2-dt4-base"
labels="dt1 dt2 dt1-base dt2-base dt1-dt3 dt2-dt4"

# running
for task in $labels; do
    {
        echo '----------------------------' Saving "$task" SurfVol Figs '----------------------------'
        fig_dir=$output_dir/SUMA_Recordings/$task
        if [[ -d "${fig_dir}" ]]; then
            rm -rf "${fig_dir}"
        fi
        mkdir -p "${fig_dir}"
        func_vol=$output_dir/$task.results/ttest+tlrc.
        afni -yesplugouts -niml $anat $func_vol &
        plugout_drive -com "QUIET_PLUGOUTS" \
            -com "SEE_OVERLAY +" \
            -com "SWITCH_UNDERLAY $anat" \
            -com "SWITCH_OVERLAY ttest 0 1" \
            -com "SET_THRESHNEW 0.001 *p" \
            -quit
        sleep 1
        # 手动点击Clusterize、NN level=2、根据Rpt中的Alpha<0.05时，确定Voxels的值
        read -r -t 600 -p "等待对AFNI进行Clusterize操作,操作后请按下回车按键以继续,否则600秒后继续!"
        plugout_drive -com "QUIET_PLUGOUTS" \
            -com "PBAR_SAVEIM $output_dir/colorbar.png" \
            -com "OPEN_WINDOW axialimage mont=4x4:10" \
            -com "OPEN_WINDOW sagittalimage mont=4x4:10" \
            -com "OPEN_WINDOW coronalimage mont=4x4:10" \
            -com "SET_DICOM_XYZ -25 40 30" \
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
        3dcalc -a $func_vol[1] -expr "ifelse(step(a*a-3.2905*3.2905),a,0)" -prefix "$fig_dir"/"$task".act.nii.gz -overwrite
    }
done

# Wait for all background tasks to complete before proceeding to the next step
wait
if [[ -L ~/temp_ln_subj_savefig ]]; then
    rm -rf ~/temp_ln_subj_savefig
fi
endTime=$(date +%s)
echo "Finished in" $((endTime - beginTime)) "s" [$begin_time - $(date)]
