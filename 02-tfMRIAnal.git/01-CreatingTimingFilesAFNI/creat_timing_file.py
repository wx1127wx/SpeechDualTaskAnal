#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2022. Institute of Health and Medical Technology, Hefei Institutes of Physical Science, CAS
# @Time     : 2022/12/20 16:28:51
# @Author   : ZL.Z
# @Email    : zzl1124@mail.ustc.edu.cn
# @Reference: None
# @FileName : creat_timing_file.py
# @Software : Python3.6; PyCharm; Windows10 / Ubuntu 18.04.5 LTS (GNU/Linux 5.4.0-79-generic x86_64)
# @Hardware : Intel Core i7-4712MQ; NVIDIA GeForce 840M / 2*X640-G30(XEON 6258R 2.7G); 3*NVIDIA GeForce RTX3090
# @Version  : V1.0 - ZL.Z：2022/12/20
# 		      First version.
# @License  : None
# @Brief    : 根据DDK双任务预实验fMRI范式所得到的的行为学数据生成AFNI格式的timing文件

import os
import shutil
import pandas as pd

if __name__ == "__main__":
    behavior_data_path = r'F:\Graduate\NeurocognitiveAssessment\认知与声音\言语与认知老化实验\言语与认知老化\data\raw_data\t-fMRI\Behavior'
    tim_res_path = r"./timing_file"
    std_timing = {"task0": ["45.62:3.2 82.2:3.2 118.82:3.2 143.22:3.2 204.23:3.2 216.43:3.2",
                            "21.22:3.2 94.42:3.2 204.22:3.2 253.03:3.2 289.63:3.2 301.83:3.2",
                            "9.03:3.2 33.42:3.2 45.62:3.2 70.02:3.2 314.02:3.2 362.82:3.2"],
                  "task1": ["21.22:3.2 94.42:3.2 106.62:3.2 167.62:3.2 265.22:3.2 338.42:3.2",
                            "45.62:3.2 106.62:3.2 216.43:3.2 277.43:3.2 326.22:3.2 338.42:3.2",
                            "82.22:3.2 143.22:3.2 167.62:3.2 179.82:3.2 253.03:3.2 301.82:3.2"],
                  "task2": ["33.42:3.2 70.02:3.2 228.63:3.2 253.03:3.2 301.82:3.2 326.22:3.2",
                            "70.02:3.2 131.02:3.2 143.22:3.2 155.42:3.2 179.82:3.2 362.82:3.2",
                            "106.62:3.2 131.02:3.2 155.42:3.2 204.23:3.2 216.43:3.2 289.63:3.2"],
                  "task3": ["179.82:3.2 192.03:3.2 240.83:3.2 289.63:3.2 350.62:3.2 362.82:3.2",
                            "82.22:3.2 167.62:3.2 192.02:3.2 228.63:3.2 265.23:3.2 314.02:3.2",
                            "57.82:3.2 94.42:3.2 228.63:3.2 265.23:3.2 277.43:3.2 326.22:3.2"],
                  "task4": ["9.03:3.2 57.82:3.2 131.02:3.2 155.42:3.2 277.43:3.2 314.02:3.2",
                            "9.03:3.2 33.43:3.2 57.82:3.2 118.82:3.2 240.83:3.2 350.62:3.2",
                            "21.23:3.2 118.82:3.2 192.02:3.2 240.83:3.2 338.42:3.2 350.62:3.2"]}
    # shutil.rmtree(tim_res_path)
    for i_subj in os.listdir(behavior_data_path):
        for t_type in ['dmBLOCK', 'BLOCK']:
            timing_save_p = os.path.join(tim_res_path, t_type, i_subj)
            if not os.path.exists(timing_save_p):
                os.makedirs(timing_save_p)
                for k_task in range(5):
                    for j_run in range(1, 4):
                        bev_f = os.path.join(behavior_data_path, i_subj, f"session_1/run{j_run}/{i_subj.split('-')[-1]}.csv")
                        if os.path.exists(bev_f):
                            bev_data = pd.read_csv(bev_f)
                            bev_data = bev_data[bev_data['trial.taskname'] == f"task{k_task}"]
                            bev_data_onset = bev_data["fix.stopped"]
                            if t_type == 'dmBLOCK':
                                bev_data_dur = bev_data["resp.started"] - bev_data_onset
                                timing_str = " ".join([f"{i}:{j}" for (i, j) in zip(bev_data_onset.tolist(),
                                                                                    bev_data_dur.tolist())])
                            else:
                                timing_str = " ".join([f"{i}" for i in bev_data_onset.tolist()])
                        else:
                            if t_type == 'dmBLOCK':
                                timing_str = std_timing[f"task{k_task}"][j_run - 1]
                            else:
                                timing_str = std_timing[f"task{k_task}"][j_run - 1].replace(":3.2", "")
                        with open(os.path.join(timing_save_p, f"timing_{['base', 'dt1', 'dt2', 'dt3', 'dt4'][k_task]}.txt"),
                                  'a', encoding='utf-8', newline='') as f:
                            f.write(timing_str + "\n")
