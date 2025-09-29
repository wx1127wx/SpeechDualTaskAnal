#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2023. Institute of Health and Medical Technology, Hefei Institutes of Physical Science, CAS
# @Time     : 2023/5/23
# @Author   : ZL.Z
# @Email    : zzl1124@mail.ustc.edu.cn
# @Reference: None
# @FileName : main.py
# @Software : Python3.6; PyCharm; Windows10
# @Hardware : Intel Core i7-4712MQ; NVIDIA GeForce 840M
# @Version  : V1.0 - ZL.Z：2023/5/23
#             First version.
# @License  : None
# @Brief    : 根据指定目录下的AFNI的预处理及GLM：afni_proc.py后的被试结果文件，自动生成3dttest++命令行

import os


if __name__ == "__main__":
    data_path = r"/home/medicaldata2/ZZLData/Datasets/audio/言语与认知老化实验/言语与认知老化/data/raw_data/t-fMRI/MRI/afni_proc_res/"
    cmd_a = " -setA $data_label_a "
    for i_dp in os.listdir(data_path):
        cmd_a += f"\\\n {i_dp.rstrip('.results')} \"$dirA/{i_dp}/stats.{i_dp.rstrip('.results')}${{use_reml}}+tlrc[$data_label_a#0_Coef]\"   "
    cmd = cmd_a
    print(cmd)
    with open("./cmd.txt", 'w') as f:
        f.write(cmd)
