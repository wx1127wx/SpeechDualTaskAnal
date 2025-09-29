# -*- coding: utf-8 -*-
# Copyright (c) 2024. Institute of Health and Medical Technology, Hefei Institutes of Physical Science, CAS
# @Time     : 2024/1/29 09:10
# @Author   : ZL.Z
# @Email    : zzl1124@mail.ustc.edu.cn
# @Reference: None
# @FileName : main.py
# @Software : Python3.6/3.9; PyCharm; Windows10 / Ubuntu 18.04.5 LTS (GNU/Linux 5.4.0-79-generic x86_64)
# @Hardware : Intel Core i7-4712MQ; NVIDIA GeForce 840M / 2*X640-G30(XEON 6258R 2.7G); 6*NVIDIA GeForce RTX3090
# @Version  : V1.0.1: 2024/10/9 小论文
#             V1.0: 2024/1/29 - 2024/1/30
#             First version.
# @License  : None
# @Brief    : DDK-WM双任务范式语音行为学数据分析

from src.config import *
import re
import datetime
import scipy.stats as stats
import scipy.signal as signal
import ptitprince as pt
from matplotlib import pyplot as plt
from matplotlib.patches import PathPatch
from matplotlib.lines import Line2D
import seaborn as sns
import pingouin as pg
from itertools import chain, combinations
from statannotations.Annotator import Annotator
from pathos.pools import ProcessPool as Pool
import math
import pywt
import glob
from collections import OrderedDict
from typing import Union, Optional, Tuple


def ttest(data_x, data_y, parametric=True, alternative='two-sided', paired=False) -> pd.DataFrame:
    """
    T检验和对应的非参数检验：单样本（data_y为一单值时）、独立样本（当数据不符合方差齐性假设时，这里自动使用Welch-t检验校正）、配对样本
    :param data_x: 输入数据，array_like
    :param data_y: 输入数据，array_like or float，当为单值时，使用单样本T检验
    :param parametric: 是否使用非参数检验, False时使用非参数检验
    :param alternative: 备择假设：双尾/单尾（'two-sided'/'greater'/'less'）
    :param paired: 是否为配对样本, True时使用配对样本检验
    :return: pandas.DataFrame，对应的检验统计结果
    """
    if parametric:  # 数据符合正态分布，参数检验：T检验
        # T检验：单样本（data_y为一单值时）、独立样本（当数据不符合方差齐性假设时，这里自动使用Welch-t检验校正）、配对样本
        res = pg.ttest(data_x, data_y, paired, alternative)
    else:  # 数据不符合正态分布，非参数检验
        if paired:  # 配对样本，使用Wilcoxon signed-rank test Wilcoxon符号秩检验
            res = pg.wilcoxon(data_x, data_y, alternative)
        else:  # 独立样本，使用Wilcoxon rank-sum test Wilcoxon秩和检验，即Mann-Whitney U Test
            res = pg.mwu(data_x, data_y, alternative)
    return res


def anova_oneway_rm(data: pd.DataFrame, save_dir: str = '', dvs=None, within=None, subject: str = '', parametric=True,
                    padjust=None, fig=True, grp_name=None, plot_kind: str = 'raincloud') -> pd.DataFrame:
    """
    单因素重复测量方差分析
    :param data: 数据，包含标签，即label
    :param save_dir: 保存结果路径
    :param dvs: 包含因变量的列的名称，列表类型
    :param within: 包含被试内组别/因子/水平的列的名称，即标签名，string类型
    :param subject: 被试标识的列的名称，string类型
    :param parametric: 是否使用非参数检验, False时使用非参数检验
    :param padjust: 事后检验的p值校正方法. None/'none': no correction; 'bonf': one-step Bonferroni correction;
                   'sidak': one-step Sidak correction; 'holm': step-down method using Bonferroni adjustments;
                   'fdr_bh': Benjamini/Hochberg FDR correction; 'fdr_by': Benjamini/Yekutieli FDR correction
    :param fig: 对否可视化各个特征的ANOVA
    :param grp_name: list或OrderedDict类型，元素为string类型,可视化中，纵坐标，即组的名称映射，顺序为从上到下的显示顺序
                list: ['0', '1', '2']，按照所列的顺序从上到下依次显示（元素必须存在于data[between]列中）
                OrderedDict: {'0':'first', '1':'second', '2':'third'}，
                将键替换为值，并按照所列的顺序从上到下依次显示（键必须存在于data[within]列中），替换后的值同时在save_dir文件中对应更改
    :param plot_kind: 绘图类型，可选violin/raincloud
    :return: 方差分析结果，包含事后检验（仅针对参数检验中的显著结果事后检验有意义）
    """
    data[within] = data[within].astype(str)
    if isinstance(grp_name, list):
        order = grp_name
    elif isinstance(grp_name, OrderedDict):
        data[within] = data[within].map(grp_name)
        order = list(grp_name.values())
    elif grp_name is None:
        order = None
    else:
        raise TypeError(f'{grp_name} 错误格式，应为list或OrderedDict格式')
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    res = pd.DataFrame()
    for dv in dvs:
        if parametric:
            # 步骤一：单因素重复测量方差分析(具体结果已包括球形检验及相应的校正)
            aov = pg.rm_anova(data, dv=dv, within=within, subject=subject, correction=True,
                              detailed=True, effsize='n2')
            # 增加效应量Effect size等级
            # J.Cohen 提出的标准:
            # 0.01- < 0.06时为小效应(small effect)，0.06 – < 0.14时为中等效应(moderate effect)，>= 0.14为大效应(large effect).
            es = aov.loc[0, 'n2']
            aov = pd.concat([aov, pd.DataFrame({'es-magnitude': ['small' if es < 0.06 else 'large'
                                                                 if es >= 0.14 else 'moderate', np.nan]})], axis=1)
        else:  # 非参数检验使用Friedman test
            # 步骤一：单因素重复测量方差分析对应的非参数检验Friedman
            aov = pg.friedman(data, dv=dv, within=within, subject=subject, method='f')
            es = aov.loc['Friedman', 'W']
            # 增加效应量Effect size等级：Friedman检验的Kendall’s W可作为效应量
            # ref: 1.	Tomczak M, Tomczak E. The need to report effect size estimates revisited.
            # An overview of some recommended measures of effect size. Trends in Sport Sciences 2014;1(21):19-25.
            # Cafiso 提出的标准 (Cafiso, S., Di Graziano, A., & Pappalardo, G. (2013). Using the Delphi method to
            # evaluate opinions of public transport managers on bus safety. Safety Science, 57, 254–263.)
            # ref: https://peterstatistics.com/CrashCourse/5-ThreeVarPair/ordinal/MultipleOrdinal3c.html
            # W <= 0.3 - Weak agreement; 0.3 < W <= 0.5 - Moderate agreement
            # 0.5 < W <= 0.7 - Good agreement; 0.7 < W - Strong agreement.
            aov = pd.concat([aov, pd.DataFrame({'es-magnitude': ['weak' if es <= 0.3 else 'large'
                                                                 if es > 0.7 else 'good'
                                                                 if 0.5 < es <= 0.7 else 'moderate']})], axis=1)
        # 步骤二：事后检验post-hoc test：配对T检验(参数检验)或Wilcoxon signed-rank test（非参数检验）
        hoc = pg.pairwise_tests(data, dv=dv, within=within, subject=subject, parametric=parametric,
                                padjust=padjust, effsize='eta-square', return_desc=True)
        p_values = [d['p-corr' if padjust not in [None, 'none'] else 'p-unc'] for i, d in hoc.iterrows()]
        ind = {}
        for i_ind in hoc.index:
            if i_ind == 0:
                ind[i_ind] = dv
            else:
                ind[i_ind] = ''
        res_dv = pd.concat([aov.add_prefix('aov-'), hoc.add_prefix('hoc-')], axis=1).rename(index=ind)
        res = pd.concat([res, res_dv])
        if fig:
            plt.figure(figsize=(7, 6), clear=True, tight_layout=True)
            plot_args = {'data': data, 'x': within, 'y': dv, 'order': order, 'orient': 'v'}
            if plot_kind == 'violin':
                ax = sns.violinplot(width=.7, cut=2, palette="pastel", scale="area", inner=None, **plot_args)
                ax = sns.swarmplot(palette="colorblind", size=3, ax=ax, zorder=5, **plot_args)
                ax = sns.boxplot(width=.15, color="gray", zorder=10, boxprops={'facecolor': 'none', "zorder": 20},
                                 showfliers=False, ax=ax, **plot_args)
                ax.set_ylim(ax.get_ylim()[0], 1.25 * ax.get_ylim()[-1])
            elif plot_kind == 'raincloud':
                plot_args['orient'] = 'h'
                ax = pt.RainCloud(cut=2, width_viol=1., pointplot=True, point_size=4, **plot_args)
            else:
                raise ValueError(f'无效的绘图类型：{plot_kind}，请从violin/raincloud中选择!')
            pairs = [(str(d['hoc-A']), str(d['hoc-B'])) for i, d in res_dv.iterrows()]
            if plot_args['orient'] == 'h':
                plot_args['x'], plot_args['y'] = dv, within
                annotator = Annotator(ax, pairs, **plot_args)
                plt.ylabel('')
                ax.xaxis.label.set_size(18)
                plt.tick_params('y', labelsize=18)
            else:
                annotator = Annotator(ax, pairs, plot='violinplot', **plot_args)
                plt.xlabel('')
                ax.yaxis.label.set_size(18)
                plt.tick_params('x', labelsize=18)
            annotator.configure(text_format='star', loc='inside', fontsize=16)
            annotator.set_pvalues(p_values).annotate()
            for sp in plt.gca().spines:
                plt.gca().spines[sp].set_color('k')
                plt.gca().spines[sp].set_linewidth(1)
            plt.xticks(fontproperties=font_family)
            plt.yticks(fontproperties=font_family)
            plt.gca().tick_params(axis='x', labelsize=14, direction='in', color='k', length=5, width=1)
            plt.gca().tick_params(axis='y', labelsize=18, direction='in', color='k', length=5, width=1)
            fig_file = os.path.join(save_dir, f"anova1_rm_{dv.split('(')[0]}.png")
            if not os.path.exists(os.path.dirname(fig_file)):
                os.makedirs(os.path.dirname(fig_file))
            plt.savefig(fig_file, dpi=600, bbox_inches='tight', pad_inches=0.02)
            plt.savefig(fig_file.replace('.png', '.tif'), dpi=600, bbox_inches='tight', pad_inches=0.02,
                        pil_kwargs={"compression": "tiff_lzw"})
            plt.show()
            plt.close()
    res.to_csv(os.path.join(save_dir, "anova1_rm.csv"), encoding="utf-8-sig")
    return res


def adjust_box_widths(ax, fac):
    """
    Adjust the widths of a seaborn-generated boxplot.
    """
    for c in ax.get_children():
        # searching for PathPatches
        if isinstance(c, PathPatch):
            # getting current width of box:
            p = c.get_path()
            verts = p.vertices
            verts_sub = verts[:-1]
            xmin = np.min(verts_sub[:, 0])
            xmax = np.max(verts_sub[:, 0])
            xmid = 0.5*(xmin+xmax)
            xhalf = 0.5*(xmax - xmin)
            # setting new width of box
            xmin_new = xmid-fac*xhalf
            xmax_new = xmid+fac*xhalf
            verts_sub[verts_sub[:, 0] == xmin, 0] = xmin_new
            verts_sub[verts_sub[:, 0] == xmax, 0] = xmax_new
            # setting new width of median line
            for l in ax.lines:
                if np.all(l.get_xdata() == [xmin, xmax]):
                    l.set_xdata([xmin_new, xmax_new])
        if isinstance(c, Line2D):
            if c.get_label() == 'cap':
                xmin = float(c.get_xdata()[0])
                xmax = float(c.get_xdata()[-1])
                xmid = 0.5 * (xmin + xmax)
                xhalf = 0.5 * (xmax - xmin)
                xmin_new = xmid - fac * xhalf
                xmax_new = xmid + fac * xhalf
                c.set_xdata([xmin_new, xmax_new])


def oculomotor_exporter(csv_speech_dir: Union[str, os.PathLike], csv_out_dir: Union[str, os.PathLike],
                        ext_all_task_from_raw: bool = False, **kwargs):
    """
    眼动记录文件和音频记录文件中提取导出对应DDK和DDKWM任务的瞳孔数据（3D的检测结果，其对测试过程中眼动仪的滑动鲁棒）
    :param csv_speech_dir: 言语测试csv文件路径
    :param csv_out_dir: 结果csv文件保存路径
    :param ext_all_task_from_raw: 是否从眼动记录原始文件（数据非常大，每个被试一个文件）中导出所有任务的瞳孔和眨眼数据（数据依然较大，每个被试一个文件）
    :param kwargs: pupil_data_exporter.PupilDataExporter参数、
                   指定从眼动记录原始文件中提取数据的被试ID列表:ext_raw_list，默认为None，即不指定，全部处理
                   眨眼指标的最小置信度阈值:min_conf_thrs_blink（默认0.2，假设人类眨眼时间是瞬时完成的，这时眨眼时长大约限制在1s左右内；
                   Pupil Core的眨眼检测算法利用一个事实:在眨眼过程中，双眼的2D瞳孔检测置信迅速下降，之后迅速上升。）
                   ext_all_task_from_raw/ext_used_task_from_exporter时并行处理的CPU核数:n_jobs
    :return: None
    """
    ddk_num, ddkwm_num, dt_num = 3, 10, 2  # 目前有3个DDK、10个DDKWM、两个双任务范式
    if kwargs.get('min_conf_thrs') is None:
        kwargs['min_conf_thrs'] = [0.8]  # 2D瞳孔检测的最小置信度阈值
    if kwargs.get('min_conf_thrs_blink') is None:
        kwargs['min_conf_thrs_blink'] = 0.2
    if kwargs.get('ext_raw_list') is None:
        kwargs['ext_raw_list'] = None
    if kwargs.get('n_jobs') is None:
        kwargs['n_jobs'] = None
    if kwargs.get('overwrite') is None:
        kwargs['overwrite'] = True
    min_conf_thrs_blink = kwargs.pop('min_conf_thrs_blink')
    ext_raw_list = kwargs.pop('ext_raw_list')
    n_jobs = kwargs.pop('n_jobs')
    overwrite = kwargs.pop('overwrite')
    assert (n_jobs is None) or (type(n_jobs) is int and n_jobs > 0) or (n_jobs == -1), 'n_jobs仅接受-1/正整数/None类型输入'
    if n_jobs == -1:
        n_jobs = None
    if ext_all_task_from_raw:
        from pupil_data_exporter import PupilDataExporter

        def ext_raw(sub):
            if (ext_raw_list is None) or (sub.split('_')[-1] in ext_raw_list):
                pup_p = os.path.join(PUPIL_RAW_PATH, f"{sub}/session_1/oculomotor_{sub.split('_')[-1]}")
                exp_p = os.path.join(PUPIL_EXPORTER_PATH, f"{sub}/session_1/oculomotor_{sub.split('_')[-1]}")
                if not os.path.exists(exp_p):
                    os.makedirs(exp_p)
                    PupilDataExporter(**kwargs).process_recordings(pup_p, csv_out_dir=exp_p,
                                                                   overwrite=overwrite, n_jobs=1)
        subj_l = [i for i in os.listdir(PUPIL_RAW_PATH) if os.path.isdir(os.path.join(PUPIL_RAW_PATH, i))]
        subj_l = sorted(subj_l)
        if n_jobs == 1:
            for i_sub in subj_l:
                ext_raw(i_sub)
        else:
            with Pool(n_jobs) as pool:
                pool.map(ext_raw, subj_l)

    def ext_exp(sub):
        print("---------- Processing %d / %d: %s ----------" % (subj_l.index(sub) + 1, len(subj_l), sub))
        exp_p = os.path.join(PUPIL_EXPORTER_PATH, f"{sub}/session_1/oculomotor_{sub.split('_')[-1]}")
        try:
            if len(os.listdir(exp_p)) > 1:  # 当一个被试的瞳孔记录文件超过一个（中断记录、多次记录），将csv文件进行拼接
                pupil_data, blink_data = pd.DataFrame(), pd.DataFrame()
                for i_rec in os.listdir(exp_p):
                    pupil_data = pd.concat([pupil_data, pd.read_csv(os.path.join(exp_p, i_rec, 'pupil.csv'))],
                                           ignore_index=True)
                    blink_data = pd.concat([blink_data, pd.read_csv(os.path.join(exp_p, i_rec, 'blinks.csv'))],
                                           ignore_index=True)
            else:
                pupil_data = pd.read_csv(os.path.join(exp_p, os.listdir(exp_p)[-1], 'pupil.csv'))
                blink_data = pd.read_csv(os.path.join(exp_p, os.listdir(exp_p)[-1], 'blinks.csv'))
        except FileNotFoundError as e:
            print(f"文件缺失/格式有误，略过：{e}")
            return None
        subj_info = pd.read_csv(os.path.join(csv_speech_dir, f"HC/{sub}/session_1/{sub.split('_')[-1]}.csv"))
        subj_info_use = subj_info.filter(regex=r'(^ddk|^ddkwm\d_)\d+\.(start|end)Time$')
        pupil_data_use, blink_data_use = pd.DataFrame(), pd.DataFrame()
        for i_ddk in range(1, ddk_num + 1):
            if pd.isnull(subj_info_use[f'ddk{i_ddk}.endTime'][0]):
                continue
            p_d = pupil_data[(pupil_data['pupil_timestamp'] >
                              float(subj_info_use[f'ddk{i_ddk}.startTime'][0].split('(')[0])-0.025) &
                             (pupil_data['pupil_timestamp'] <
                              float(subj_info_use[f'ddk{i_ddk}.endTime'][0].split('(')[0])+0.02) &
                             (pupil_data['model_confidence'] == 1.0) & (pupil_data['method'] == 'pye3d 0.3.0 real-time')]
            p_d_use = p_d[['pupil_timestamp', 'eye_id', 'diameter',
                           'diameter_3d', 'model_confidence', 'confidence']].reset_index(drop=True)
            i_pupil_data_use = pd.concat([pd.DataFrame({'id': [sub.split('_')[-1]]*p_d_use.index.__len__(),
                                                        'task': [f'ddk{i_ddk}']*p_d_use.index.__len__()}),
                                          p_d_use], axis=1)
            pupil_data_use = pd.concat([pupil_data_use, i_pupil_data_use], ignore_index=True)
            b_d = blink_data[(blink_data['start_timestamp'] >
                              float(subj_info_use[f'ddk{i_ddk}.startTime'][0].split('(')[0]) - 0.025) &
                             (blink_data['start_timestamp'] <
                              float(subj_info_use[f'ddk{i_ddk}.endTime'][0].split('(')[0]) + 0.02) &
                             (blink_data['confidence'] >= min_conf_thrs_blink)]
            b_d_use = b_d[['start_timestamp', 'duration_s', 'confidence']].reset_index(drop=True)
            i_blink_data_use = pd.concat([pd.DataFrame({'id': [sub.split('_')[-1]] * b_d_use.index.__len__(),
                                                        'task': [f'ddk{i_ddk}'] * b_d_use.index.__len__()}),
                                          b_d_use], axis=1)
            blink_data_use = pd.concat([blink_data_use, i_blink_data_use], ignore_index=True)
        for i_dt in range(1, dt_num + 1):
            for i_ddkwm in range(1, ddkwm_num + 1):
                if pd.isnull(subj_info_use[f'ddkwm{i_dt}_{i_ddkwm}.endTime'][0]):
                    continue
                p_d = pupil_data[(pupil_data['pupil_timestamp'] >
                                  float(subj_info_use[f'ddkwm{i_dt}_{i_ddkwm}.startTime'][0].split('(')[0])-0.025) &
                                 (pupil_data['pupil_timestamp'] <
                                  float(subj_info_use[f'ddkwm{i_dt}_{i_ddkwm}.endTime'][0].split('(')[0])+0.02) &
                                 (pupil_data['model_confidence'] == 1.0) & (pupil_data['method'] == 'pye3d 0.3.0 real-time')]
                p_d_use = p_d[['pupil_timestamp', 'eye_id', 'diameter',
                               'diameter_3d', 'model_confidence', 'confidence']].reset_index(drop=True)
                i_pupil_data_use = pd.concat([pd.DataFrame({'id': [sub.split('_')[-1]]*p_d_use.index.__len__(),
                                                            'task': [f'ddkwm{i_dt}_{i_ddkwm}']*p_d_use.index.__len__()}),
                                              p_d_use], axis=1)
                pupil_data_use = pd.concat([pupil_data_use, i_pupil_data_use], ignore_index=True)
                b_d = blink_data[(blink_data['start_timestamp'] >
                                  float(subj_info_use[f'ddkwm{i_dt}_{i_ddkwm}.startTime'][0].split('(')[0]) - 0.025) &
                                 (blink_data['start_timestamp'] <
                                  float(subj_info_use[f'ddkwm{i_dt}_{i_ddkwm}.endTime'][0].split('(')[0]) + 0.02) &
                                 (blink_data['confidence'] >= min_conf_thrs_blink)]
                b_d_use = b_d[['start_timestamp', 'duration_s', 'confidence']].reset_index(drop=True)
                i_blink_data_use = pd.concat([pd.DataFrame({'id': [sub.split('_')[-1]] * b_d_use.index.__len__(),
                                                            'task': [f'ddkwm{i_dt}_{i_ddkwm}'] * b_d_use.index.__len__()}),
                                              b_d_use], axis=1)
                blink_data_use = pd.concat([blink_data_use, i_blink_data_use], ignore_index=True)
        return pupil_data_use, blink_data_use
    subj_l = [i for i in os.listdir(PUPIL_EXPORTER_PATH) if os.path.isdir(os.path.join(PUPIL_EXPORTER_PATH, i))]
    subj_l = sorted(subj_l)
    if n_jobs == 1:
        res = []
        for i_sub in subj_l:
            res.append(ext_exp(i_sub))
    else:
        with Pool(n_jobs) as pool:
            res = pool.map(ext_exp, subj_l)
    pupil_data_use_all, blink_data_use_all = pd.DataFrame(), pd.DataFrame()
    for _res in res:
        if _res is not None:
            pupil_data_use_all = pd.concat([pupil_data_use_all, _res[0]], ignore_index=True)
            blink_data_use_all = pd.concat([blink_data_use_all, _res[-1]], ignore_index=True)
    if pupil_data_use_all.empty:
        print("无瞳孔数据！")
        return None
    pupil_data_use_all.sort_values(by=['id', 'pupil_timestamp', 'eye_id'], inplace=True)
    pupil_data_use_all.to_csv(os.path.join(csv_out_dir, 'pupil_data_use_raw.csv'), encoding="utf-8-sig", index=False)
    blink_data_use_all.sort_values(by=['id', 'start_timestamp'], inplace=True)
    blink_data_use_all.to_csv(os.path.join(csv_out_dir, 'blink_data_use_raw.csv'), encoding="utf-8-sig", index=False)


def pupil_preprocessing(pupil_data_raw: pd.DataFrame, csv_speech_dir: Union[str, os.PathLike],
                        csv_out_dir: Union[str, os.PathLike], run_standard_prep: bool = False,
                        n_jobs=None) -> pd.DataFrame:
    """
    全部被试瞳孔数据标准化预处理及排除缺失数据过多的被试或试次
    ref: M. E. Kret and E. E. Sjak-Shie, "Preprocessing pupil size data: Guidelines and code,"
         Behavior Research Methods, vol. 51, no. 3, pp. 1336-1342, 2018.
         https://dr-jt.github.io/pupillometry/articles/preprocess_overview.html
    :param pupil_data_raw: 全部被试原始瞳孔数据
    :param csv_speech_dir: 言语测试csv文件路径
    :param csv_out_dir: 保存的预处理结果csv文件路径
    :param run_standard_prep: 是否从所导出的所有任务的瞳孔原始数据csv文件中预处理数据，否则直接导入已处理好的pupil_data_prep.csv
    :param n_jobs: 并行核数
    :return: pd.DataFrame，预处理及排除缺失数据过多的被试或试次后的瞳孔结果
    """
    subjs = pupil_data_raw['id'].unique().tolist()
    eyes = pupil_data_raw['eye_id'].unique().tolist()
    tasks = pupil_data_raw['task'].unique().tolist()
    if run_standard_prep:

        def pupil_prep(sub_id):
            print("---------- Processing %d / %d: %s ----------" % (subjs.index(sub_id) + 1, len(subjs), sub_id))
            pd_pupil_pre_sub = pd.DataFrame()
            for eye_id in eyes:
                for task in tasks:
                    # Step1: 删除眨眼附近的瞳孔数据(该步在Pupil导出数据时通过瞳孔置信度已排除)和
                    # 排除伪迹(一方面Pupil Labs使用的3D瞳孔物理模型一定程度上避免了伪迹；另一方面，使用绝对中位值偏差MAD排除残余伪迹)
                    pd_data_use_all = pupil_data_raw[(pupil_data_raw['id'] == sub_id)
                                                     & (pupil_data_raw['eye_id'] == eye_id)
                                                     & (pupil_data_raw['task'] == task)]
                    pd_data = pd_data_use_all['diameter_3d']
                    if len(pd_data) > 1:
                        mad = stats.median_abs_deviation(pd_data)
                        pd_data_no_art = pd_data.where((pd_data - pd_data.median()).abs() <= 16 * mad)
                        pd_no_art = pd.concat([pd_data_use_all.drop(columns=['diameter_3d', 'model_confidence']),
                                               pd_data_no_art], axis=1)
                        # Step2: 线性插值以重采样替换nan值（根据参考文献，连续缺失值超过250ms时不做插值）
                        pd_no_art['diameter_3d'] = pd_no_art['diameter_3d'].interpolate(limit=int(0.25 * 120))
                        # Step3: 应用零相位数字低通滤波器进行平滑，截止频率为4Hz
                        ba = signal.butter(N=3, Wn=4, fs=120, btype='low', analog=False, output='ba')
                        try:
                            pd_no_art['diameter_3d'] = signal.filtfilt(ba[0], ba[1], pd_no_art['diameter_3d'])
                            pd_pupil_pre = pd_no_art.copy()
                        except ValueError:
                            pd_pupil_pre = pd_no_art.copy()
                    else:
                        pd_pupil_pre = pd_data_use_all.drop(columns=['model_confidence'])
                    pd_pupil_pre_sub = pd.concat([pd_pupil_pre_sub, pd_pupil_pre])
            subj_info = pd.read_csv(glob.glob(os.path.join(csv_speech_dir, f"HC/*_{sub_id}/session_1/{sub_id}.csv"))[0])
            subj_info_fil = subj_info.filter(regex=r'number|(^ddk|^ddkwm\d_)\d+\.(start|end)Time$')
            subj_info_use = subj_info_fil.astype({'number': int})
            subj_info_use.set_index('number', inplace=True)
            events = []
            np_filter_all = pd_pupil_pre_sub.to_numpy()
            for i_rec in range(len(pd_pupil_pre_sub)):
                rec_task = np_filter_all[i_rec, 1]
                rec_ts = np_filter_all[i_rec, 2]
                if rec_task.startswith('ddkwm'):
                    t_s = float(
                        subj_info_use.at[sub_id, f"ddkwm{rec_task[5]}_{rec_task[7:]}.startTime"].split('(')[0]) - 0.025
                    t_e = float(subj_info_use.at[sub_id, f"ddkwm{rec_task[5]}_{rec_task[7:]}.endTime"].split('(')[0]) + 0.02
                    if t_s <= rec_ts < t_s + 0.6:
                        event = 'prepare'
                    elif t_s + 0.6 <= rec_ts < t_s + 1.2:
                        event = '+'
                    elif t_s + 1.2 <= rec_ts < t_s + 1.4:
                        event = 'beforeblank'
                    elif (t_s + 2.0 <= rec_ts < t_s + 2.2) or (t_s + 2.8 <= rec_ts < t_s + 3.0) or \
                            (t_s + 3.6 <= rec_ts < t_s + 3.8):
                        event = 'isi'
                    elif (t_s + 1.4 <= rec_ts < t_s + 2.0) or (t_s + 2.2 <= rec_ts < t_s + 2.8) or \
                            (t_s + 3.0 <= rec_ts < t_s + 3.6) or (t_s + 3.8 <= rec_ts < t_s + 4.4):
                        event = 'encode'
                    elif t_s + 4.4 <= rec_ts < t_e:
                        event = 'afterblank'
                    else:
                        event = 'none'
                else:
                    event = 'ddk'
                events.append(event)
            pd_pupil_pre_sub.insert(2, 'event', events)
            pd_pupil_pre_sub.sort_values(by=['id', 'pupil_timestamp', 'eye_id'], inplace=True, ignore_index=True)
            return pd_pupil_pre_sub

        assert (n_jobs is None) or (type(n_jobs) is int and n_jobs > 0) or (
                n_jobs == -1), 'n_jobs仅接受-1/正整数/None类型输入'
        if n_jobs == -1:
            n_jobs = None
        if n_jobs == 1:
            res = []
            for i_subj in subjs:
                res.append(pupil_prep(i_subj))
        else:
            with Pool(n_jobs) as pool:
                res = pool.map(pupil_prep, subjs)
        pd_pupil_pre_all = pd.DataFrame()
        for _res in res:
            pd_pupil_pre_all = pd.concat([pd_pupil_pre_all, _res])
        pd_pupil_pre_all.sort_values(by=['id', 'pupil_timestamp', 'eye_id'], inplace=True, ignore_index=True)
        pd_pupil_pre_all.to_csv(os.path.join(csv_out_dir, 'pupil_data_prep.csv'), encoding="utf-8-sig", index=False)
    else:
        pd_pupil_pre_all = pd.read_csv(os.path.join(csv_out_dir, 'pupil_data_prep.csv'))
    # 排除无效被试：缺失值过多或瞳孔数据过少的被试（标准差法）
    p_id_grp = pd_pupil_pre_all.dropna(subset=['diameter_3d']).groupby(['id'], as_index=False).count()
    p_num = p_id_grp['diameter_3d']
    num_thr_low = p_num.mean() - 2 * p_num.std()
    id_filter = p_id_grp[p_num >= num_thr_low]['id']
    pupil_filter_subj = pd_pupil_pre_all[pd_pupil_pre_all['id'].isin(id_filter)]
    # 排除被试中的无效试次/事件：被试中的缺失值过多或瞳孔数据过少的试次（标准差法）
    p_id_task_grp = pupil_filter_subj.dropna(subset=['diameter_3d']).groupby(['id', 'eye_id', 'task', 'event'],
                                                                             as_index=False).count()
    pupil_filter = pd.DataFrame()
    for i_eye in pupil_filter_subj['eye_id'].unique().tolist():
        for j_task in pupil_filter_subj['task'].unique().tolist():
            for k_event in pupil_filter_subj['event'].unique().tolist():
                p_num = p_id_task_grp[(p_id_task_grp['eye_id'] == i_eye) & (p_id_task_grp['task'] == j_task) &
                                      (p_id_task_grp['event'] == k_event)]['diameter_3d']
                num_thr_low = p_num.mean() - 2 * p_num.std()
                id_filter = p_id_task_grp[(p_id_task_grp['eye_id'] == i_eye) & (p_id_task_grp['task'] == j_task) &
                                          (p_id_task_grp['event'] == k_event) & (p_num >= num_thr_low)]['id']
                pupil_task = pupil_filter_subj[(pupil_filter_subj['eye_id'] == i_eye) & (pupil_filter_subj['task'] == j_task)
                                               & (pupil_filter_subj['event'] == k_event)]
                pupil_filter_task = pupil_task[pupil_task['id'].isin(id_filter)]
                pupil_filter = pd.concat([pupil_filter, pupil_filter_task])
                # print(i_eye, j_task, k_event, p_num.mean(), num_thr_low)
    pupil_filter.sort_values(by=['id', 'pupil_timestamp', 'eye_id'], inplace=True, ignore_index=True)
    pupil_filter.to_csv(os.path.join(csv_out_dir, 'pupil_data.csv'), encoding="utf-8-sig", index=False)
    return pupil_filter


def get_lhipa(pupil_data: pd.DataFrame, csv_out_dir: Union[str, os.PathLike], run_lhipa: bool = False,
              n_jobs=None) -> pd.DataFrame:
    """
    计算所需的全部被试的LHIPA
    :param pupil_data: 瞳孔数据
    :param csv_out_dir: 结果保存路径
    :param run_lhipa: 是否运行计算lhipa，否则从已计算的csv文件导入
    :param n_jobs: 计算lhipa的并行核数
    :return: pd.DataFrame，所需的全部被试的LHIPA
    """
    subjs = pupil_data['id'].unique().tolist()
    eyes = pupil_data['eye_id'].unique().tolist()
    tasks = pupil_data['task'].unique().tolist()
    if run_lhipa:

        def calu_lhipa(sub_id):
            print("---------- Processing %d / %d: %s ----------" % (subjs.index(sub_id) + 1, len(subjs), sub_id))
            pd_lhipa = pd.DataFrame()
            for eye_id in eyes:
                for task in tasks:
                    for_calu = pupil_data[(pupil_data['id'] == sub_id) & (pupil_data['eye_id'] == eye_id) &
                                          (pupil_data['task'] == task)]
                    if task.startswith('ddkwm'):
                        baseline = for_calu[for_calu['event'].isin(['prepare', '+', 'beforeblank'])]
                        stimulation = for_calu[for_calu['event'].isin(['encode', 'isi', 'afterblank'])]
                        # print(eye_id, task, len(baseline), len(stimulation))
                        # 缺失数据的阈值10-30%（ref: Winn, et al., Trends Hear, 2018; Geller, et al., Behav. Res. Methods, 2020）
                        if len(baseline) < 0.2 * int(np.ceil((0.6 + 0.6 + 0.2) * 120)):
                            lhipa_val_base = np.nan
                        else:
                            lhipa_val_base = lhipa(baseline, 120) if lhipa(baseline, 120) else np.nan
                        if len(stimulation) < 0.2 * int(np.ceil((0.6*4 + 0.2*3 + 0.6) * 120)):
                            lhipa_val_stim = np.nan
                        else:
                            lhipa_val_stim = lhipa(stimulation, 120) if lhipa(stimulation, 120) else np.nan
                        pd_lhipa = pd.concat([pd_lhipa,
                                              pd.DataFrame({'id': [sub_id]*2, 'eye_id': [eye_id]*2, 'task': [task]*2,
                                                            'event': ['baseline', 'stimulation'],
                                                            'lhipa': [lhipa_val_base, lhipa_val_stim]})])
                    else:
                        if len(for_calu) < 0.2 * int(np.ceil(5.0 * 120)):
                            lhipa_val_ddk = np.nan
                        else:
                            lhipa_val_ddk = lhipa(for_calu, 120) if lhipa(for_calu, 120) else np.nan
                        pd_lhipa = pd.concat([pd_lhipa, pd.DataFrame({'id': [sub_id], 'eye_id': [eye_id], 'task': [task],
                                                                      'event': ['ddk'], 'lhipa': [lhipa_val_ddk]})])
            return pd_lhipa
        assert (n_jobs is None) or (type(n_jobs) is int and n_jobs > 0) or (
                n_jobs == -1), 'n_jobs仅接受-1/正整数/None类型输入'
        if n_jobs == -1:
            n_jobs = None
        if n_jobs == 1:
            res = []
            for i_subj in subjs:
                res.append(calu_lhipa(i_subj))
        else:
            with Pool(n_jobs) as pool:
                res = pool.map(calu_lhipa, subjs)
        pd_lhipa_all = pd.DataFrame()
        for _res in res:
            pd_lhipa_all = pd.concat([pd_lhipa_all, _res])
        pd_lhipa_all['task'] = pd_lhipa_all['task'].astype('category').cat. \
            set_categories([f'ddk{i}' for i in range(1, 4)] + [f'ddkwm1_{i}' for i in range(1, 11)] +
                           [f'ddkwm2_{i}' for i in range(1, 11)])
        pd_lhipa_all['event'] = pd_lhipa_all['event'].astype('category').cat.set_categories(['ddk', 'baseline', 'stimulation'])
        pd_lhipa_all.sort_values(by=['id', 'eye_id', 'task', 'event'], inplace=True, ignore_index=True)
        pd_lhipa_all.to_csv(os.path.join(csv_out_dir, 'lhipa_all.csv'), encoding="utf-8-sig", index=False)
    else:
        pd_lhipa_all = pd.read_csv(os.path.join(csv_out_dir, 'lhipa_all.csv'))
    # print(pd_lhipa_all)
    pd_lhipa_nonan = pd_lhipa_all[~pd_lhipa_all['task'].str.contains(r'ddk\d')].dropna(subset=['lhipa']).reset_index(drop=True)
    pupil_lhipa = pd.DataFrame()
    task_map = {r'ddkwm1_\d+': 'ddkwm1', r'ddkwm2_\d+': 'ddkwm2'}
    for i_sub in pd_lhipa_nonan['id'].unique().tolist():
        lhipa_sub = pd_lhipa_nonan[pd_lhipa_nonan['id'] == i_sub]
        lhipa_z = stats.zscore(lhipa_sub['lhipa'])
        lhipa_sub.insert(len(lhipa_sub.columns), 'lhipa_z', lhipa_z)
        for i_task in [r'ddkwm1_\d+', r'ddkwm2_\d+']:
            for i_event in ['baseline', 'stimulation']:
                lhipa_sub_task = lhipa_sub[lhipa_sub['task'].str.contains(i_task) & (lhipa_sub['event'] == i_event)]
                if not lhipa_sub_task.empty:
                    lhipa_sub_task_filter = lhipa_sub_task.groupby(['eye_id'], as_index=False).mean()
                    lhipa_sub_task_eye0 = lhipa_sub_task_filter[lhipa_sub_task_filter['eye_id'] == 0]
                    lhipa_sub_task_eye1 = lhipa_sub_task_filter[lhipa_sub_task_filter['eye_id'] == 1]
                    lhipa_eye0 = lhipa_sub_task_eye0['lhipa'].values[-1] if len(lhipa_sub_task_eye0['lhipa'].values) else np.nan
                    lhipa_eye1 = lhipa_sub_task_eye1['lhipa'].values[-1] if len(lhipa_sub_task_eye1['lhipa'].values) else np.nan
                    lhipa_eye_mean = np.nanmean([lhipa_eye0, lhipa_eye1])
                    lhipa_eye_diff = abs(lhipa_eye0 - lhipa_eye1) if (not np.isnan(lhipa_eye0)) and \
                                                                     (not np.isnan(lhipa_eye1)) else np.nan
                    lhipa_z_eye0 = lhipa_sub_task_eye0['lhipa_z'].values[-1] if len(lhipa_sub_task_eye0['lhipa_z'].values) else np.nan
                    lhipa_z_eye1 = lhipa_sub_task_eye1['lhipa_z'].values[-1] if len(lhipa_sub_task_eye1['lhipa_z'].values) else np.nan
                    lhipa_z_eye_mean = np.nanmean([lhipa_z_eye0, lhipa_z_eye1])
                    lhipa_z_eye_diff = abs(lhipa_z_eye0 - lhipa_z_eye1) if (not np.isnan(lhipa_z_eye0)) and \
                                                                           (not np.isnan(lhipa_z_eye1)) else np.nan
                    pupil_lhipa = pd.concat([pupil_lhipa,
                                             pd.DataFrame({'id': [i_sub], 'task': [task_map[i_task]], 'event': [i_event],
                                                           'lhipa_eye0': [lhipa_eye0],
                                                           'lhipa_eye1': [lhipa_eye1], 'lhipa_eye_mean': [lhipa_eye_mean],
                                                           'lhipa_eye_diff': [lhipa_eye_diff],
                                                           'lhipa_z_eye0': [lhipa_z_eye0], 'lhipa_z_eye1': [lhipa_z_eye1],
                                                           'lhipa_z_eye_mean': [lhipa_z_eye_mean],
                                                           'lhipa_z_eye_diff': [lhipa_z_eye_diff]})])
                else:
                    pupil_lhipa = pd.concat([pupil_lhipa,
                                             pd.DataFrame({'id': [i_sub], 'task': [task_map[i_task]], 'event': [i_event],
                                                           'lhipa_eye0': [np.nan], 'lhipa_eye1': [np.nan],
                                                           'lhipa_eye_mean': [np.nan], 'lhipa_eye_diff': [np.nan],
                                                           'lhipa_z_eye0': [np.nan], 'lhipa_z_eye1': [np.nan],
                                                           'lhipa_z_eye_mean': [np.nan], 'lhipa_z_eye_diff': [np.nan]})])
    pupil_lhipa['event'] = pupil_lhipa['event'].astype('category').cat.set_categories(['baseline', 'stimulation'])
    pupil_lhipa.sort_values(by=['id', 'task', 'event'], inplace=True, ignore_index=True)
    pupil_lhipa.to_csv(os.path.join(csv_out_dir, 'pupil_lhipa.csv'), encoding="utf-8-sig", index=False)
    # print(pupil_lhipa)
    return pupil_lhipa


def get_blink(blink_data: pd.DataFrame, csv_speech_dir: Union[str, os.PathLike],
              csv_out_dir: Union[str, os.PathLike]) -> pd.DataFrame:
    """
    计算所需的全部被试的眨眼率指标
    :param blink_data: 眨眼数据
    :param csv_speech_dir: 言语测试csv文件路径
    :param csv_out_dir: 结果保存路径
    :return: pd.DataFrame，所需的全部被试的眨眼率指标
    """
    pd_blink = blink_data.copy()
    subjs = pd_blink['id'].unique().tolist()
    subj_info_use = pd.DataFrame()
    for i_csv in glob.glob(os.path.join(csv_speech_dir, f"HC/*/session_1/*.csv")):
        subj_info = pd.read_csv(i_csv)
        subj_info_fil = subj_info.filter(regex=r'number|(^ddk|^ddkwm\d_)\d+\.(start|end)Time$')
        subj_info_fil = subj_info_fil.astype({'number': int})
        subj_info_use = pd.concat([subj_info_use, subj_info_fil])
    subj_info_use.set_index('number', inplace=True)
    events = []
    np_filter_all = pd_blink.to_numpy()
    for i_rec in range(len(pd_blink)):
        rec_id = np_filter_all[i_rec, 0]
        rec_task = np_filter_all[i_rec, 1]
        rec_ts = np_filter_all[i_rec, 2]
        if rec_task.startswith('ddkwm'):
            t_s = float(
                subj_info_use.at[rec_id, f"ddkwm{rec_task[5]}_{rec_task[7:]}.startTime"].split('(')[0]) - 0.025
            t_e = float(subj_info_use.at[rec_id, f"ddkwm{rec_task[5]}_{rec_task[7:]}.endTime"].split('(')[0]) + 0.02
            if t_s <= rec_ts < t_s + 0.6:
                event = 'prepare'
            elif t_s + 0.6 <= rec_ts < t_s + 1.2:
                event = '+'
            elif t_s + 1.2 <= rec_ts < t_s + 1.4:
                event = 'beforeblank'
            elif (t_s + 2.0 <= rec_ts < t_s + 2.2) or (t_s + 2.8 <= rec_ts < t_s + 3.0) or \
                    (t_s + 3.6 <= rec_ts < t_s + 3.8):
                event = 'isi'
            elif (t_s + 1.4 <= rec_ts < t_s + 2.0) or (t_s + 2.2 <= rec_ts < t_s + 2.8) or \
                    (t_s + 3.0 <= rec_ts < t_s + 3.6) or (t_s + 3.8 <= rec_ts < t_s + 4.4):
                event = 'encode'
            elif t_s + 4.4 <= rec_ts < t_e:
                event = 'afterblank'
            else:
                event = 'none'
        else:
            event = 'ddk'
        events.append(event)
    pd_blink.insert(2, 'event', events)
    pd_blink.sort_values(by=['id', 'start_timestamp'], inplace=True, ignore_index=True)
    show_seq = {'prepare': 0.6, '+': 0.6, 'beforeblank': 0.2, 'isi': 0.6, 'afterblank': 0.6, 'encode': 2.4}
    ddk_num, ddkwm_num = 3, 10  # 目前有3个DDK、10个DDKWM
    blink_event_detail = pd.DataFrame()
    for sg in show_seq.keys():
        blink_use = pd_blink[pd_blink['event'].isin(['ddk', sg])]
        for i_task in [r'ddk\d', r'ddkwm1_\d+', r'ddkwm2_\d+']:
            blink = blink_use[blink_use['task'].str.contains(i_task)]
            dur_c = blink.groupby(['id'], as_index=False).count()['task']
            # 计算眨眼率：每分钟眨眼次数（每个任务包含多次测试）
            blink_ratio = dur_c / (5.0 * ddk_num if i_task == r'ddk\d' else show_seq[sg] * ddkwm_num) * 60
            blink_ratio = pd.concat([pd.DataFrame({'id': blink['id'].unique()}), blink_ratio], axis=1)
            blink_ratio.rename(columns={'task': 'blink_ratio (cpm)'}, inplace=True)
            blink_task = pd.DataFrame({'id': subjs}).merge(blink_ratio, how='outer', on='id').reset_index(drop=True)
            blink_task.insert(1, 'task', [i_task.split('\d')[0].rstrip('_')] * len(blink_task))
            if i_task == r'ddk\d':
                blink_task.insert(2, 'event', ['ddk'] * len(blink_task))
            else:
                blink_task.insert(2, 'event', [sg] * len(blink_task))
            blink_event_detail = pd.concat([blink_event_detail, blink_task])
    blink_event_detail.drop_duplicates(inplace=True)
    blink_event_detail['event'] = blink_event_detail['event'].astype('category').cat.\
        set_categories(['ddk', 'prepare', '+', 'beforeblank', 'isi', 'encode', 'afterblank'])
    blink_event_detail.sort_values(by=['id', 'task', 'event'], inplace=True, ignore_index=True)
    blink_event_detail.to_csv(os.path.join(csv_out_dir, 'blink_event_detail.csv'), encoding="utf-8-sig", index=False)
    blink_event = blink_event_detail[~blink_event_detail['task'].str.contains(r'ddk\d')]
    baseline = blink_event[blink_event['event'].isin(['prepare', '+', 'beforeblank'])]
    stimulation = blink_event[blink_event['event'].isin(['encode', 'isi', 'afterblank'])]
    base_grp = baseline.groupby(['id', 'task'], as_index=False).mean()
    base_grp.insert(2, 'event', ['baseline'] * len(base_grp))
    stim_grp = stimulation.groupby(['id', 'task'], as_index=False).mean()
    stim_grp.insert(2, 'event', ['stimulation'] * len(stim_grp))
    pupil_blink = pd.concat([base_grp, stim_grp])
    pupil_blink['event'] = pupil_blink['event'].astype('category').cat.set_categories(['baseline', 'stimulation'])
    pupil_blink.sort_values(by=['id', 'task', 'event'], inplace=True, ignore_index=True)
    pupil_blink.to_csv(os.path.join(csv_out_dir, 'pupil_blink.csv'), encoding="utf-8-sig", index=False)
    return pupil_blink


def lhipa(pupil_data: pd.DataFrame, fs: int = None):
    """
     计算The Low/High Index of Pupillary Activity，LHIPA
     ref: A. T. Duchowski, et al., "The Index of Pupillary Activity," presented at the Proceedings of the 2018 CHI
     Conference on Human Factors in Computing Systems, 2018.
     A. T. Duchowski, et al., "The Low/High Index of Pupillary Activity," presented at the Proceedings of the 2020 CHI
     Conference on Human Factors in Computing Systems, 2020.
    :param pupil_data: 瞳孔数据
    :param fs: 采样率
    :return: LHIPA，The Low/High Index of Pupillary Activity
    """
    d = pupil_data['diameter_3d'].to_numpy()
    # get signal duration (in seconds)
    tt = pupil_data['pupil_timestamp'].to_numpy()[-1] - pupil_data['pupil_timestamp'].to_numpy()[0]
    if fs:
        # print(tt, len(d), len(pupil_data['pupil_timestamp'].to_numpy()))
        resampled_x, resampled_t = signal.resample(d, int(np.ceil(tt * fs)), pupil_data['pupil_timestamp'].to_numpy())
        d = resampled_x.copy()
        tt = resampled_t[-1] - resampled_t[0]
    # find max decomposition level
    w = pywt.Wavelet('sym4')
    maxlevel = pywt.dwt_max_level(len(d), filter_len=w.dec_len)
    # set high and low frequency band indeces, ie. j
    hif, lof = 1, int(maxlevel / 2)
    # print(hif, lof, maxlevel, len(d), w.dec_len)
    # get detail coefficients of pupil diameter signal d
    try:
        cD_H = pywt.downcoef('d', d, 'sym4', 'per', level=hif)
        cD_L = pywt.downcoef('d', d, 'sym4', 'per', level=lof)
    except ValueError:
        return np.nan
    # normalize by 1/√(2j)
    cD_H[:] = [x / math.sqrt(2 ** hif) for x in cD_H]
    cD_L[:] = [x / math.sqrt(2 ** lof) for x in cD_L]
    # print('cD_H: ', cD_H[:5])
    # print('cD_L: ', cD_L[:5])
    # obtain the LH:HF ratio
    cD_LH = cD_L.copy()
    for i in range(len(cD_L)):
        cD_LH[i] = cD_L[i] / cD_H[int(((2**lof)/(2**hif)) * i)]
    # detect modulus maxima, see Duchowski et al. [15]
    cD_LHm = modmax(cD_LH)
    # print('cD_LH: ', cD_LH[:5])
    # threshold using universal threshold λuniv = σ√(2logn), where σ is the standard deviation of the noise
    lambda_univ = np.std(cD_LHm) * math.sqrt(2.0 * np.log2(len(cD_LHm)))
    # print(lambda_univ)
    cD_LHt = pywt.threshold(cD_LHm, lambda_univ, mode="less")
    # print(cD_LHt)
    # compute LHIPA: count of coefficients per second
    ctr = (np.fabs(cD_LHt) > 0).sum()
    lhipa_val = ctr / tt
    return lhipa_val


def modmax(d):
    """
     Modulus maxima detection
     ref: A. T. Duchowski, et al., "The Index of Pupillary Activity," presented at the Proceedings of the 2018 CHI
     Conference on Human Factors in Computing Systems, 2018.
     A. T. Duchowski, et al., "The Low/High Index of Pupillary Activity," presented at the Proceedings of the 2020 CHI
     Conference on Human Factors in Computing Systems, 2020.
    :param d: 1D signal
    :return: Modulus maxima
    """
    # compute signal modulus
    m = [0.0] * len(d)
    for i in range(len(d)):
        m[i] = math.fabs(d[i])
    # if value is larger than both neighbours, and strictly larger than either, then it is a local maximum
    t = [0.0] * len(d)
    for i in range(len(d)):
        ll = m[i - 1] if i >= 1 else m[i]
        oo = m[i]
        rr = m[i+1] if i < len(d)-2 else m[i]
        if (ll <= oo and oo >= rr) and (ll < oo or oo > rr):
            # compute magnitude
            t[i] = math.sqrt(d[i]**2)
        else:
            t[i] = 0.0
    return t


def task_pupil_compare(data_acoustic: pd.DataFrame, data_oculomotor: pd.DataFrame, save_dir: Union[str, os.PathLike],
                       wm_thr: int = 5, outliers_detect: bool = True, parametric: Optional[bool] = None,
                       between_dual: bool = False, fig: bool = True,
                       plot_kind: str = 'raincloud') -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    对双任务的眼动特征进行T检验并绘图:证明双任务的负荷增加能够反映在眼动特征上
    :param data_acoustic: 全部被试的声学指标数据，至少需要包含WM-Correct Count数据
    :param data_oculomotor: 全部被试的眼动指标数据
    :param save_dir: 结果保存路径
    :param wm_thr: 工作记忆次任务的正确数阈值
                   即最少wm_thr个回答正确的被试才被作为待检验的条目（总共10个试次，不少于wm_thr个准确回答才作为有效的双任务范式数据）
    :param outliers_detect: 是否进行异常值检测（对每列特征的异常值赋值为nan）
    :param parametric: T检验/对应的非参数检验，True/False/None，None为根据每次检验的数据的正态性进行自动判断
    :param between_dual: 对否计算和显示不同难度双任务之间的差异
    :param fig: 对否可视化各个特征的t检验结果
    :param plot_kind: 绘图类型，可选violin/raincloud
    :return: t检验结果
    """
    dt_num = 2  # 目前有两个双任务范式
    data_ocu = data_oculomotor.copy()
    for i_dual in range(1, dt_num + 1):
        valid_id = data_acoustic[data_acoustic[f'13_DDKWM-T{i_dual}-WM-Correct Count'] >= wm_thr]['id']
        data_ocu[data_oculomotor.drop(columns=['id', 'task', 'event']).columns] = \
            data_ocu[data_oculomotor.drop(columns=['id', 'task', 'event']).columns].mask(~data_ocu['id'].isin(valid_id))
    data_ocu.columns = data_ocu.columns.str.replace(r' \(.*\)$', '', regex=True)
    res_dict_de, res_dict_p = {}, {}
    for i_feat in data_ocu.drop(columns=['id', 'task', 'event']).columns:
        res_task_de, res_task_p, dt_data = [], [], []
        for j_dual in range(1, dt_num + 1):
            base_sel = data_ocu[(data_ocu['task'] == f'ddkwm{j_dual}') &
                                (data_ocu['event'] == 'baseline')][['id', i_feat]].reset_index(drop=True)
            stim_sel = data_ocu[(data_ocu['task'] == f'ddkwm{j_dual}') &
                                (data_ocu['event'] == 'stimulation')][['id', i_feat]].reset_index(drop=True)
            # 利用箱线图法，即根据四分位数进行异常值检测
            if outliers_detect:
                base_thr_high = base_sel.quantile(0.75) + 1.5 * (base_sel.quantile(0.75) - base_sel.quantile(0.25))
                base_thr_low = base_sel.quantile(0.25) - 1.5 * (base_sel.quantile(0.75) - base_sel.quantile(0.25))
                stim_thr_high = stim_sel.quantile(0.75) + 1.5 * (stim_sel.quantile(0.75) - stim_sel.quantile(0.25))
                stim_thr_low = stim_sel.quantile(0.25) - 1.5 * (stim_sel.quantile(0.75) - stim_sel.quantile(0.25))
                # 异常值则判为无效数据，此时清空，即赋值为nan
                base_sel_filter = base_sel[(base_sel >= base_thr_low) & (base_sel <= base_thr_high)]
                stim_sel_filter = stim_sel[(stim_sel >= stim_thr_low) & (stim_sel <= stim_thr_high)]
            else:
                base_sel_filter = base_sel.copy()
                stim_sel_filter = stim_sel.copy()
            x_data = base_sel_filter[i_feat]
            y_data = stim_sel_filter[i_feat]
            dt_data.append(y_data)
            if parametric is None:  # 根据正态性自动判断
                stat, p = stats.shapiro(x_data - y_data)  # 正态性检验：Shapiro-Wilk Test
                if p < 0.05:  # 非正态分布，使用非参数检验
                    parametric = False
                else:
                    parametric = True
            t_res = ttest(x_data, y_data, parametric, alternative='greater', paired=True)
            par_dict = {}
            for k_par in t_res.columns:
                par_dict[k_par] = t_res[k_par][0]
            res_task_de.append(par_dict)
            res_task_p.append(t_res['p-val'][0])
        if between_dual:
            for _x_data, _y_data in combinations(dt_data, 2):
                _xy_data = pd.concat([_x_data, _y_data], axis=1, join='inner')
                x_data, y_data = _xy_data.iloc[:, 0], _xy_data.iloc[:, 1]
                if parametric is None:  # 根据正态性自动判断
                    stat, p = stats.shapiro(x_data - y_data)  # 正态性检验：Shapiro-Wilk Test
                    if p < 0.05:  # 非正态分布，使用非参数检验
                        parametric = False
                    else:
                        parametric = True
                t_res = ttest(x_data, y_data, parametric, alternative='greater', paired=True)
                par_dict = {}
                for k_par in t_res.columns:
                    par_dict[k_par] = t_res[k_par][0]
                res_task_de.append(par_dict)
                res_task_p.append(t_res['p-val'][0])
        res_dict_de[i_feat] = res_task_de
        res_dict_p[i_feat] = res_task_p
    ind = ['baseline-low', 'baseline-high']
    if between_dual:
        ind += ['low-high']
    p_value_de = pd.DataFrame(res_dict_de, index=ind)
    p_value = pd.DataFrame(res_dict_p, index=ind)
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    p_value_de.to_csv(os.path.join(save_dir, 'task_pupil_ttest_detail.csv'), encoding="utf-8-sig")
    p_value.to_csv(os.path.join(save_dir, 'task_pupil_ttest_p.csv'), encoding="utf-8-sig")
    print(p_value)
    if fig:
        data_sel = data_ocu.drop(columns=['id'])
        data_sel['task'] = data_sel['task'].str.replace(r'^ddkwm1$', 'Low Difficulty', regex=True).\
            str.replace(r'^ddkwm2$', 'High Difficulty', regex=True)
        data_sel['event'] = data_sel['event'].str.replace('stimulation', 'Load Stage', regex=True).\
            str.replace('baseline', 'Baseline Stage', regex=True)
        for i_col in data_sel.drop(columns=['task', 'event']).columns:
            lab = 'LHIPA' if i_col.startswith('lhipa') else 'Blink Rate' if i_col.startswith('blink') else i_col
            plt.figure(figsize=(7, 6), clear=True, tight_layout=True)
            plot_args = {'data': data_sel, 'x': 'task', 'y': i_col, 'hue': 'event',
                         'order': ['Low Difficulty', 'High Difficulty'], 'hue_order': ['Baseline Stage', 'Load Stage'], 'orient': 'v'}
            if plot_kind == 'violin':
                ax = sns.violinplot(width=.7, cut=2, palette="pastel", scale="area", inner=None, **plot_args)
                ax = sns.boxplot(width=.7, color="gray", zorder=10, boxprops={'facecolor': 'none', "zorder": 20},
                                 showfliers=True, ax=ax, capprops={'label': 'cap'}, **plot_args)
                ax = pt.stripplot(width=.7, size=3, dodge=True, ax=ax, zorder=5, **plot_args)
                adjust_box_widths(ax, .2)
                ax.set_ylim(ax.get_ylim()[0], 1.25 * ax.get_ylim()[-1])
                lg_loc = 'upper left'
            elif plot_kind == 'raincloud':
                plot_args['orient'] = 'h'
                ax = pt.RainCloud(cut=2, width_viol=1., pointplot=True, point_size=4,
                                  alpha=.65, dodge=True, move=.2, **plot_args)
                lg_loc = 'upper right'
            else:
                raise ValueError(f'无效的绘图类型：{plot_kind}，请从violin/raincloud中选择!')
            ax.xaxis.label.set_size(18)
            ax.yaxis.label.set_size(18)
            plt.xticks(fontproperties=font_family)
            plt.yticks(fontproperties=font_family)
            plt.gca().tick_params(labelsize=14, direction='in', color='k', length=5, width=1)
            pairs = [(('Low Difficulty', 'Baseline Stage'), ('Low Difficulty', 'Load Stage')),
                     (('High Difficulty', 'Baseline Stage'), ('High Difficulty', 'Load Stage'))]
            if between_dual:
                pairs += [(('Low Difficulty', 'Load Stage'), ('High Difficulty', 'Load Stage'))]
            pv = p_value[i_col]
            if plot_args['orient'] == 'h':
                plot_args['x'], plot_args['y'] = i_col, 'task'
                annotator = Annotator(ax, pairs, **plot_args)
                plt.xlabel(lab)
                plt.ylabel('')
                ax.xaxis.label.set_size(18)
                plt.tick_params('y', labelsize=18)
            else:
                annotator = Annotator(ax, pairs, plot='violinplot', **plot_args)
                plt.xlabel('')
                plt.ylabel(lab)
                ax.yaxis.label.set_size(18)
                plt.tick_params('x', labelsize=18)
            annotator.configure(text_format='star', loc='inside', fontsize=18)
            annotator.set_pvalues(pv).annotate()
            _handles, _labels = ax.get_legend_handles_labels()
            handles, labels = [], []
            for i_it in range(len(_labels)):
                if not isinstance(_handles[i_it], Line2D):
                    handles.append(_handles[i_it])
                if _labels[i_it] != 'cap':
                    labels.append(_labels[i_it])
            plt.legend(handles=handles[:len(set(labels))], labels=labels[:len(set(labels))], loc=lg_loc,
                       prop={'family': font_family, 'size': 14}, title=None)
            for sp in plt.gca().spines:
                plt.gca().spines[sp].set_color('k')
                plt.gca().spines[sp].set_linewidth(1)
            fig_file = os.path.join(save_dir, f"{i_col}.png")
            if not os.path.exists(os.path.dirname(fig_file)):
                os.makedirs(os.path.dirname(fig_file))
            plt.savefig(fig_file, dpi=600, bbox_inches='tight', pad_inches=0.02)
            plt.savefig(fig_file.replace('.png', '.tif'), dpi=600, bbox_inches='tight', pad_inches=0.02,
                        pil_kwargs={"compression": "tiff_lzw"})
            plt.show()
            plt.close()
    return p_value, p_value_de


def task_acoustic_compare(data_acoustic: pd.DataFrame, save_dir: Union[str, os.PathLike], wm_thr: int = 5,
                          outliers_detect: bool = True, parametric: Optional[bool] = None,
                          padjust: Optional[str] = None, fig: bool = True, grp_name=None,
                          plot_kind: str = 'raincloud') -> pd.DataFrame:
    """
    对3个不同任务的声学特征进行单因素重复测量方差分析并绘图:证明双任务的负荷增加能够反映在声学特征上
    :param data_acoustic: 全部数据
    :param save_dir: 结果保存路径
    :param wm_thr: 工作记忆次任务的正确数阈值
                   即最少wm_thr个回答正确的被试才被作为待检验的条目（总共10个试次，不少于wm_thr个准确回答才作为有效的双任务范式数据）
    :param outliers_detect: 是否进行异常值检测（对每列特征的异常值赋值为nan）
    :param parametric: 单因素重复测量方差分析对应的非参数检验(Friedman test)，True/False/None(根据每个水平数据的正态性自动判断)
    :param padjust: 事后检验的p值校正方法
    :param fig: 对否可视化统计检验结果
    :param grp_name: list或OrderedDict类型，元素为string类型,可视化中，纵坐标，即组的名称映射，顺序为从上到下的显示顺序
                list: ['0', '1', '2']，按照所列的顺序从上到下依次显示（元素必须存在于data[between]列中）
                OrderedDict: {'0':'first', '1':'second', '2':'third'}，
                将键替换为值，并按照所列的顺序从上到下依次显示（键必须存在于data[within]列中），替换后的值同时在save_dir文件中对应更改
    :param plot_kind: 绘图类型，可选violin/raincloud
    :return: 统计检验结果
    """
    dt_num = 2  # 目前有两个双任务范式
    pd_task = pd.DataFrame(list(chain(*[['Sig'] * data_acoustic.shape[0]] +
                                       [[f'Dual{i}'] * data_acoustic.shape[0] for i in range(1, dt_num + 1)])),
                           columns=['task'])
    pd_data = pd.concat([data_acoustic.filter(regex=r'id|^02_DDK-').rename(columns=lambda x: re.sub(r'^02_DDK-', '', x))]
                        + [data_acoustic.filter(regex=f'id|^13_DDKWM-T{i}-DDK-').rename(columns=lambda x:
                        re.sub(f'^13_DDKWM-T{i}-DDK-', '', x)) for i in range(1, dt_num + 1)]).reset_index(drop=True)
    pd_data.insert(1, 'task', pd_task['task'])
    pd_data['task'] = pd_data['task'].astype('category').cat.set_categories(['Sig'] + [f'Dual{i}' for i in range(1, dt_num + 1)])
    pd_data.sort_values(by=['id', 'task'], ignore_index=True, inplace=True)
    pd_data.columns = pd_data.columns.str.replace(r' \(.*\)$', '', regex=True).str.replace('pause', 'Pause').\
        str.replace('rate', 'Rate').str.replace('duration', 'Duration').str.replace('regularity', 'Regularity')
    data_aco = pd_data.copy()
    for i_dual in range(1, dt_num + 1):
        valid_id = data_acoustic[data_acoustic[f'13_DDKWM-T{i_dual}-WM-Correct Count'] >= wm_thr]['id']
        data_aco[pd_data.drop(columns=['id', 'task']).columns] = data_aco[pd_data.drop(columns=['id', 'task']).columns].\
            mask(~data_aco['id'].isin(valid_id))
    data_aco_filter = data_aco.copy()
    for i_aco in data_aco_filter.drop(columns=['id', 'task']).columns:
        # 利用箱线图法，即根据四分位数进行异常值检测：Q1-1.5*IQR（下四分位数减去1.5倍四分位距）
        if outliers_detect:
            data_aco_grp = data_aco_filter.groupby(['task'])[i_aco]
            thr_high = data_aco_grp.quantile(0.75) + 1.5 * (data_aco_grp.quantile(0.75) - data_aco_grp.quantile(0.25))
            thr_low = data_aco_grp.quantile(0.25) - 1.5 * (data_aco_grp.quantile(0.75) - data_aco_grp.quantile(0.25))
            for j_task in thr_high.index:  # 异常值则判为无效数据，此时清空，即赋值为nan
                data_aco_filter[i_aco] = data_aco_filter[i_aco].mask((data_aco_filter['task'] == j_task) &
                                                                     ((data_aco_filter[i_aco] < thr_low[j_task]) |
                                                                     (data_aco_filter[i_aco] > thr_high[j_task])))
        if parametric is None:  # 根据正态性自动判断
            norm = True  # 正态性检验：Shapiro-Wilk Test
            for j_task in data_aco_filter['task'].unique().tolist():
                stat, p = stats.shapiro(data_aco_filter[data_aco_filter['task'] == j_task][i_aco])
                if p < 0.05:  # 非正态分布
                    norm = False
                    break
            if norm:
                parametric = True
            else:  # 非正态分布，使用非参数检验
                parametric = False
    aov = anova_oneway_rm(data=data_aco_filter, dvs=data_aco_filter.drop(columns=['id', 'task']).columns,
                          within='task', subject='id', parametric=parametric, save_dir=save_dir, padjust=padjust,
                          fig=fig, plot_kind=plot_kind, grp_name=grp_name)
    return aov


def wm_compare(data_acoustic: pd.DataFrame, save_dir: Union[str, os.PathLike], outliers_detect: bool = True,
               parametric: Optional[bool] = None, fig: bool = True) -> pd.DataFrame:
    """
    对双任务的行为学指标（正确数）进行配对样本T检验并绘图
    :param data_acoustic: 全部被试的声学指标数据
    :param save_dir: 结果保存路径
    :param outliers_detect: 是否进行异常值检测（对每列特征的异常值赋值为nan）
    :param parametric: 是否使用非参数检验, False时使用非参数检验（Wilcoxon signed-rank test Wilcoxon符号秩检验）
    :param fig: 对否可视化统计检验结果
    :return: 配对样本T检验结果
    """
    dt_num = 2  # 目前有两个双任务范式
    cognitive_load = pd.DataFrame()
    for i_dual in range(1, dt_num + 1):
        aco = data_acoustic.iloc[:, data_acoustic.columns.str.contains(f'id|13_DDKWM-T{i_dual}-WM-Correct Count')]
        aco.columns = aco.columns.str.replace(f'13_DDKWM-T{i_dual}-', '').str.replace(r' \(.*\)$', '', regex=True).\
            str.replace('pause', 'Pause').str.replace('rate', 'Rate').str.replace('duration', 'Duration').\
            str.replace('regularity', 'Regularity')
        aco.insert(1, 'task', [f'ddkwm{i_dual}'] * len(aco))
        cognitive_load = pd.concat([cognitive_load, aco])
    cognitive_load.sort_values(by=['id', 'task'], ignore_index=True, inplace=True)
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    cognitive_load.to_csv(os.path.join(save_dir, f'cognitive_load_dtc_acoustic.csv'), encoding="utf-8-sig", index=False)
    res_tt = pd.DataFrame()
    for i_feat in cognitive_load.drop(columns=['id', 'task']).columns:
        cl_filter = cognitive_load.copy()
        if outliers_detect:  # 利用箱线图法，即根据四分位数进行异常值检测
            data_grp = cl_filter.groupby(['task'])[i_feat]
            thr_high = data_grp.quantile(0.75) + 1.5 * (data_grp.quantile(0.75) - data_grp.quantile(0.25))
            thr_low = data_grp.quantile(0.25) - 1.5 * (data_grp.quantile(0.75) - data_grp.quantile(0.25))
            for j_task in thr_high.index:  # 异常值则判为无效数据，此时清空，即赋值为nan
                cl_filter[i_feat] = cl_filter[i_feat].mask((cl_filter['task'] == j_task) &
                                                           ((cl_filter[i_feat] < thr_low[j_task]) |
                                                            (cl_filter[i_feat] > thr_high[j_task])))
        x_data = cl_filter[cl_filter['task'] == 'ddkwm1'][i_feat]
        y_data = cl_filter[cl_filter['task'] == 'ddkwm2'][i_feat]
        if parametric is None:  # 根据正态性自动判断
            stat, p = stats.shapiro(x_data - y_data)  # 正态性检验：Shapiro-Wilk Test
            if p < 0.05:  # 非正态分布，使用非参数检验
                parametric = False
            else:
                parametric = True
        alternative = {'WM-Correct Count': 'greater'}
        t_res = ttest(x_data, y_data, parametric, alternative=alternative[i_feat], paired=True)
        res_tt = pd.concat([res_tt, t_res.rename(index={'T-test': f'Low-High({i_feat})'})])
        res_tt.to_csv(os.path.join(save_dir, "dtc_acoustic_ttest.csv"), encoding="utf-8-sig")
        if fig:
            plt.figure(figsize=(7, 6), clear=True, tight_layout=True)
            annotator = Annotator(ax=plt.gca(), pairs=[('ddkwm1', 'ddkwm2')], data=cl_filter, x='task', y=i_feat)
            ax = pg.plot_paired(data=cl_filter, dv=i_feat, within='task', subject='id', boxplot_in_front=True, ax=plt.gca(),
                                boxplot_kwargs={'palette': "pastel", 'width': 0.2}, pointplot_kwargs={'scale': 0.3})
            ax.set_xticklabels(['Low Difficulty', 'High Difficulty'])
            plt.xlabel('')
            ax.yaxis.label.set_size(18)
            plt.tick_params('x', labelsize=18)
            annotator.configure(text_format='star', loc='inside', fontsize=16)
            annotator.set_pvalues(t_res['p-val']).annotate()
            for sp in plt.gca().spines:
                plt.gca().spines[sp].set_color('k')
                plt.gca().spines[sp].set_linewidth(1)
            plt.xticks(fontproperties=font_family)
            plt.yticks(fontproperties=font_family)
            plt.gca().tick_params(axis='x', labelsize=18, direction='in', color='k', length=5, width=1)
            plt.gca().tick_params(axis='y', labelsize=14, direction='in', color='k', length=5, width=1)
            fig_file = os.path.join(save_dir, f"dtc_ttest_{i_feat}.png")
            if not os.path.exists(os.path.dirname(fig_file)):
                os.makedirs(os.path.dirname(fig_file))
            plt.savefig(fig_file, dpi=600, bbox_inches='tight', pad_inches=0.02)
            plt.savefig(fig_file.replace('.png', '.tif'), dpi=600, bbox_inches='tight', pad_inches=0.02,
                        pil_kwargs={"compression": "tiff_lzw"})
            plt.show()
            plt.close()
    print(res_tt)
    return res_tt


if __name__ == "__main__":
    start_time = datetime.datetime.now()
    print(f"---------- Start Time ({os.path.basename(__file__)}): {start_time.strftime('%Y-%m-%d %H:%M:%S')} ----------")
    current_path = os.path.dirname(os.path.realpath(__file__))
    data_path = os.path.join(current_path, r"data")
    res_path = os.path.join(current_path, r"results")
    subinfo_speech_f = os.path.join(data_path, 'subinfo_speech.csv')

    run_ocu_exp_flag = False  # STEP1: 眼动瞳孔数据导出
    if run_ocu_exp_flag:
        oculomotor_exporter(PREP_PATH, csv_out_dir=data_path, ext_all_task_from_raw=False, n_jobs=-1, overwrite=False)
    run_ocu_ext_flag = False  # STEP2: 眼动特征提取
    if run_ocu_ext_flag:
        # 瞳孔数据预处理
        pupil_prep_filter = pupil_preprocessing(pd.read_csv(os.path.join(data_path, 'pupil_data_use_raw.csv')),
                                                PREP_PATH, data_path, run_standard_prep=False, n_jobs=-1)
        # 计算瞳孔数据的认知负荷指标：LHIPA
        data_lhipa = get_lhipa(pupil_prep_filter, data_path, run_lhipa=False, n_jobs=-1)
        # 计算眨眼数据的认知负荷指标：眨眼率
        data_blink = get_blink(pd.read_csv(os.path.join(data_path, 'blink_data_use_raw.csv')), PREP_PATH, data_path)
    else:
        data_lhipa = pd.read_csv(os.path.join(data_path, 'pupil_lhipa.csv'))
        data_blink = pd.read_csv(os.path.join(data_path, 'pupil_blink.csv'))
    subinfo_speech = pd.read_csv(subinfo_speech_f)
    acoustic_feat = subinfo_speech[(subinfo_speech['label'] == 'aging') & (subinfo_speech['session'] == 1)]
    acoustic_feat = acoustic_feat.dropna(subset=acoustic_feat.filter(regex="^(02_|13_)").columns)
    less_than_30 = acoustic_feat[acoustic_feat['age'] < 30].sample(n=30, random_state=rs)
    acoustic_feat_match = pd.concat([acoustic_feat[acoustic_feat['age'] >= 30], less_than_30]).sort_index(ignore_index=True)
    # print(acoustic_feat_match.shape, acoustic_feat_match[acoustic_feat_match['age'] > 50].shape)  # 101, 52
    # acoustic_feat_match.to_csv(os.path.join(data_path, "data_match30.csv"), encoding="utf-8-sig", index=False)
    cl_pupil = data_lhipa[['id', 'task', 'event', 'lhipa_z_eye_mean']].merge(data_blink, how='outer',
                                                                             on=['id', 'task', 'event'])
    cl_pupil.sort_values(by=['id', 'task', 'event'], inplace=True, ignore_index=True)
    run_ocu_comp_flag = False  # STEP3: 双任务测量认知负荷的有效性验证（不同难度任务的瞳孔数据认知负荷指标的差异性比较）
    if run_ocu_comp_flag:
        task_pupil_compare(acoustic_feat, cl_pupil, os.path.join(res_path, 'task_pupil_compare'),
                           wm_thr=5, outliers_detect=True,
                           parametric=True, between_dual=True, fig=True, plot_kind='violin')
        task_pupil_compare(acoustic_feat_match, cl_pupil, os.path.join(res_path, 'matched/task_pupil_compare'),
                           wm_thr=5, outliers_detect=True,
                           parametric=True, between_dual=True, fig=True, plot_kind='violin')
    run_aco_comp_flag = False  # STEP4: 声学双任务指标衡量认知负荷的有效性验证（双任务与基线的声学指标差异比较）
    if run_aco_comp_flag:
        # 双任务与基线的声学指标差异性比较(python3.9)
        task_acoustic_compare(acoustic_feat, os.path.join(res_path, 'task_acoustic_compare'), wm_thr=5,
                              outliers_detect=True, parametric=True, padjust='fdr_bh', fig=True, plot_kind='raincloud',
                              grp_name=OrderedDict({'Sig': 'Single-Task', 'Dual1': 'Low Difficulty',
                                                    'Dual2': 'High Difficulty'}))
        # 不同难度次任务WM准确数的差异性比较
        wm_compare(acoustic_feat, os.path.join(res_path, 'wm_compare'), outliers_detect=True, parametric=True, fig=True)

        task_acoustic_compare(acoustic_feat_match, os.path.join(res_path, 'matched/task_acoustic_compare'), wm_thr=5,
                              outliers_detect=True, parametric=True, padjust='fdr_bh', fig=True, plot_kind='raincloud',
                              grp_name=OrderedDict({'Sig': 'Single-Task', 'Dual1': 'Low Difficulty',
                                                    'Dual2': 'High Difficulty'}))
        wm_compare(acoustic_feat_match, os.path.join(res_path, 'matched/wm_compare'),
                   outliers_detect=True, parametric=True, fig=True)

    end_time = datetime.datetime.now()
    print(f"---------- End Time ({os.path.basename(__file__)}): {end_time.strftime('%Y-%m-%d %H:%M:%S')} ----------")
    print(f"---------- Time Used ({os.path.basename(__file__)}): {end_time - start_time} ----------")
    with open(os.path.join(current_path, "results/finished.txt"), "w") as ff:
        ff.write(f"------------------ Started at {start_time.strftime('%Y-%m-%d %H:%M:%S')} "
                 f"({os.path.basename(__file__)}) -------------------\r\n")
        ff.write(f"------------------ Finished at {end_time.strftime('%Y-%m-%d %H:%M:%S')} "
                 f"({os.path.basename(__file__)}) -------------------\r\n")
        ff.write(f"------------------ Time Used {end_time - start_time} "
                 f"({os.path.basename(__file__)}) -------------------\r\n")
