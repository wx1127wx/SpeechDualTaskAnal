# -*- coding: utf-8 -*-
# Copyright (c) 2021. Institute of Health and Medical Technology, Hefei Institutes of Physical Science, CAS
# @Time     : 2021/8/13 15:17 
# @Author   : ZL.Z
# @Email    : zzl1124@mail.ustc.edu.cn
# @Reference: None
# @FileName : config.py
# @Software : Python3.6; PyCharm; Windows10
# @Hardware : Intel Core i7-4712MQ; NVIDIA GeForce 840M
# @Version  : V1.0 - ZL.Z：2021/8/13 - 2021/8/13
# 		      First version.
# @License  : None
# @Brief    : 配置文件

import os
import random
import platform
import pandas as pd
import numpy as np
import matplotlib
import warnings

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 200)
pd.set_option('display.width', 200)
np.set_printoptions(suppress=True, formatter={'all': lambda x: str(x)})
os.environ['OUTDATED_IGNORE'] = '1'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
# warnings.filterwarnings('ignore', message=r'Support for FigureCanvases without a required_interactive_framework.')

if platform.system() == 'Windows':
    PREP_PATH = r"F:\Graduate\NeurocognitiveAssessment\认知与声音\言语与认知老化实验\言语与认知老化\analysis\preprocessed_data"
    PUPIL_RAW_PATH = r'F:\Graduate\NeurocognitiveAssessment\认知与声音\言语与认知老化实验\言语与认知老化\data\raw_data\pupil\raw_data'
    PUPIL_EXPORTER_PATH = r'F:\Graduate\NeurocognitiveAssessment\认知与声音\言语与认知老化实验\言语与认知老化\data\raw_data\pupil\exporter_data'
else:
    PREP_PATH = r"/home/zlzhang/data/言语与认知老化实验/言语与认知老化/analysis/preprocessed_data"
    PUPIL_RAW_PATH = r'/home/medicaldata/ZZLData/Datasets/audio/言语与认知老化实验/言语与认知老化/data/raw_data/pupil/raw_data'
    PUPIL_EXPORTER_PATH = r'/home/medicaldata/ZZLData/Datasets/audio/言语与认知老化实验/言语与认知老化/data/raw_data/pupil/exporter_data'
font_family = 'Times New Roman'
matplotlib.rcParams["font.family"] = font_family


def setup_seed(seed: int):
    """
    全局固定随机种子
    :param seed: 随机种子值
    :return: None
    """
    random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)
    np.random.seed(seed)


rs = 323
setup_seed(rs)
