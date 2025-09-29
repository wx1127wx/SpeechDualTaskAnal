# -*- coding: utf-8 -*-
# Copyright (c) 2024. Institute of Health and Medical Technology, Hefei Institutes of Physical Science, CAS
# @Time     : 2024/2/5 11:46
# @Author   : ZL.Z
# @Email    : zzl1124@mail.ustc.edu.cn
# @Reference: None
# @FileName : main.py
# @Software : Python3.6; PyCharm / VS Code; Windows10 / Ubuntu 18.04.5 LTS (GNU/Linux 5.4.0-79-generic x86_64)
# @Hardware : Intel Core i7-4712MQ; NVIDIA GeForce 840M / 2*X640-G30(XEON 6258R 2.7G); 6*NVIDIA GeForce RTX3090
# @Version  : V1.1.0: 2024/11/27 线性回归改为SVM，并用AIC指标
#             V1.0.1: 2024/10/9 小论文
#             V1.0: 2024/2/5 - 2024/2/18
#             First version.
# @License  : None
# @Brief    : 言语双任务评估认知老化

import os
import random
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch
from matplotlib.lines import Line2D
import matplotlib.ticker as ticker
import ptitprince as pt
from statannotations.Annotator import Annotator
import seaborn as sns
import scikit_posthocs as scp
import datetime
from typing import Union, Optional, Tuple
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import mean_squared_error
from scipy.stats import chi2_contingency
from itertools import combinations
import warnings
os.environ['OUTDATED_IGNORE'] = '1'
import pingouin as pg

pd.set_option('display.max_columns', 15)
pd.set_option('display.max_rows', 10)
pd.set_option('display.width', 200)
np.set_printoptions(suppress=True, formatter={'all': lambda x: str(x)})
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
font_family = 'Times New Roman'
matplotlib.rcParams["font.family"] = font_family
warnings.filterwarnings("ignore", message="Invalid x-position found.*", category=UserWarning)


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


def get_demographic_info_task_status(data_dem_speech: pd.DataFrame, save_dir: Union[str, os.PathLike],
                                     min_age: Optional[int] = None, age_interval: int = 10, parametric: bool = True,
                                     equ_var: bool = True, padjust: Optional[str] = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    获取各年龄段人口统计学信息，并进行组间比较；以及获取各任务完成情况
    :param data_dem_speech: 包含人口统计学信息及言语特征的数据
    :param save_dir: 保存结果路径
    :param min_age: 年龄段的最小其实年龄
    :param age_interval: 每个年龄区间段的大小
    :param parametric: 是否使用非参数检验, False时使用非参数检验
    :param equ_var: 当使用参数检验时, 是否根据方差齐性检验来选择事后检验的方法
    :param padjust: 使用非参数检验时事后检验的p值校正方法. None/'none': no correction; 'bonf': one-step Bonferroni correction;
                   'sidak': one-step Sidak correction; 'holm': step-down method using Bonferroni adjustments;
                   'fdr_bh': Benjamini/Hochberg FDR correction; 'fdr_by': Benjamini/Yekutieli FDR correction
    :return: 各年龄段人口统计学信息、各任务完成情况
    """
    # STEP1: 获取各年龄段人口统计学信息，并进行组间比较
    min_age = data_dem_speech['age'].min() if min_age is None else min_age
    max_age = data_dem_speech['age'].max()
    bins = list(range(min_age, max_age + age_interval, age_interval))
    labels = [f'{start} ~ {start + age_interval - 1}' for start in bins[:-1]]
    labels[-1] = f'≥ {bins[-2]}'
    dem_used = data_dem_speech.copy()[['age', 'sex', 'edu', 'MMSE', 'DST', 'IQCODE']]
    dem_used['age_group'] = pd.cut(dem_used['age'], bins=bins, labels=labels, right=False)
    _dem_used = dem_used.copy()
    _dem_used['age_group'] = '全部样本'
    dem_used = pd.concat([dem_used, _dem_used])
    pd_dem = pd.DataFrame({'项目名': ['样本数 (男/女)', '年龄 (年)', '教育水平 (年)', 'MMSE (分)', 'DST (分)', 'IQCODE (分)']})
    grp_sex, grp_age, grp_edu, grp_mmse, grp_dst, grp_iqcode = pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), \
        pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    for grp_name, grp_data in dem_used.groupby('age_group'):
        if grp_name != '全部样本':
            grp_sex = pd.concat([grp_sex, pd.DataFrame({grp_name: [grp_data[grp_data['sex'] == 0]['sex'].count(), 
                                                                   grp_data[grp_data['sex'] == 1]['sex'].count()]})],
                                axis=1)
            grp_age = pd.concat([grp_age, pd.DataFrame({grp_name: grp_data['age'].to_list()})], axis=1)
            grp_edu = pd.concat([grp_edu, pd.DataFrame({grp_name: grp_data['edu'].to_list()})], axis=1)
            grp_mmse = pd.concat([grp_mmse, pd.DataFrame({grp_name: grp_data['MMSE'].to_list()})], axis=1)
            grp_dst = pd.concat([grp_dst, pd.DataFrame({grp_name: grp_data['DST'].to_list()})], axis=1)
            grp_iqcode = pd.concat([grp_iqcode, pd.DataFrame({grp_name: grp_data['IQCODE'].to_list()})], axis=1)
        samp_num = f"{grp_data.shape[0]} ({grp_data[grp_data['sex'] == 1]['sex'].count()}" \
            f" / {grp_data[grp_data['sex'] == 0]['sex'].count()})"
        age = f"{grp_data['age'].mean():.2f} ± {grp_data['age'].std():.2f}" \
            f" ({grp_data['age'].min()} ~ {grp_data['age'].max()})"
        edu = f"{grp_data['edu'].mean():.2f} ± {grp_data['edu'].std():.2f}" \
            f" ({grp_data['edu'].min():.2f} ~ {grp_data['edu'].max():.2f})"
        mmse = f"{grp_data['MMSE'].mean():.2f} ± {grp_data['MMSE'].std():.2f}" \
            f" ({grp_data['MMSE'].min()} ~ {grp_data['MMSE'].max()})"
        dst = f"{grp_data['DST'].mean():.2f} ± {grp_data['DST'].std():.2f}" \
            f" ({grp_data['DST'].min()} ~ {grp_data['DST'].max()})"
        iqcode = f"{grp_data['IQCODE'].mean():.2f} ± {grp_data['IQCODE'].std():.2f}" \
            f" ({grp_data['IQCODE'].min():.2f} ~ {grp_data['IQCODE'].max():.2f})"
        pd_dem = pd.concat([pd_dem, pd.DataFrame({grp_name: [samp_num, age, edu, mmse, dst, iqcode]})], axis=1)
    pd_dem.insert(1, '全部样本', pd_dem.pop('全部样本'))
    
    def _anova_oneway(data: pd.DataFrame):
        dv, between = 'val', 'grp'
        data_long = data.melt(var_name=between, value_name=dv)
        if parametric:
            # 步骤一：方差齐性检验
            data_nonan = data_long.dropna(subset=[dv])  # pg.homoscedasticity要求非nan
            homo = pg.homoscedasticity(data_nonan, dv=dv, group=between)  # Levene方法
            # 满足或设定不用方差齐性，使用标准的ANOVA和Tukey-HSD post-hoc test
            if (not equ_var) or (equ_var and homo['equal_var'].tolist()[0]):
                # 步骤二：方差分析(ANOVA)
                aov = pg.anova(data_long, dv=dv, between=between)
                # 步骤三：事后检验post-hoc test：Tukey-HSD post-hoc test
                hoc = pg.pairwise_tukey(data_long, dv=dv, between=between, effsize='eta-square')
                p_corrected = 'p-tukey'
            else:  # 不满足方差齐性，使用Welch ANOVA和Games-Howell post-hoc test
                # 步骤二：方差分析(ANOVA)
                aov = pg.welch_anova(data_long, dv=dv, between=between)
                # 步骤三：事后检验post-hoc test：Games-Howell post-hoc test
                hoc = pg.pairwise_gameshowell(data_long, dv=dv, between=between, effsize='eta-square')
                p_corrected = 'pval'
        else:  # 非参数检验直接使用Kruskal-Wallis秩和检验
            # 步骤一：方差分析(ANOVA)
            aov = pg.kruskal(data_long, dv=dv, between=between).rename(columns={'H': 'F'})
            # 步骤二：事后检验post-hoc test：Conover´s post-hoc test
            _data_long = data_long.dropna(subset=[dv])  # scikit_posthocs包需提前去除np.nan才能保持结果与R或SPSS一致
            hoc_p = scp.posthoc_conover(_data_long, dv, between, p_adjust=None).to_numpy()
            hoc_p_adj = scp.posthoc_conover(_data_long, dv, between, p_adjust=padjust).to_numpy()
            pval = hoc_p[np.triu_indices_from(hoc_p, k=1)]
            pval_adj = hoc_p_adj[np.triu_indices_from(hoc_p_adj, k=1)]
            _grp = _data_long.groupby(between, observed=True)[dv]
            _labels = np.array(list(_grp.groups.keys()))
            gmeans = _grp.mean().to_numpy()
            ng = aov['ddof1'][0] + 1
            g1, g2 = np.array(list(combinations(np.arange(ng), 2))).T
            hoc = pd.DataFrame({'A': _labels[g1], 'B': _labels[g2], 'mean(A)': gmeans[g1],
                                'mean(B)': gmeans[g2], 'pval': pval, 'pval_adj': pval_adj})
            p_corrected = 'pval_adj' if padjust is not None else 'pval'
        _f_val, _p_val = aov['F'][0], aov['p-unc'][0]
        _hoc_p_val = {}
        for i_comp in hoc.index:
            _hoc_p_val[f"{hoc.loc[i_comp, 'A']} vs. {hoc.loc[i_comp, 'B']}"] = hoc.loc[i_comp, p_corrected]
        return _f_val, _p_val, _hoc_p_val
    
    chi, p_sex, dof, exp = chi2_contingency(grp_sex, correction=False)
    comp_grp, comp_hoc = [f"χ2 = {chi:.3f}, p = {p_sex:.3f}"], ["—"]
    for item in grp_age, grp_edu, grp_mmse, grp_dst, grp_iqcode:
        f_val, p_val, hoc_p_val = _anova_oneway(item)
        comp_grp.append(f"F = {f_val:.3f}, p = {p_val:.3f}")
        hoc_t = ''
        for key, val in hoc_p_val.items():
            hoc_t += f"{key}: {val:.3f}, "
        comp_hoc.append(hoc_t)
    pd_dem = pd.concat([pd_dem, pd.DataFrame({'组比较': comp_grp, '事后比较': comp_hoc})], axis=1)
    pd_dem.to_csv(os.path.join(save_dir, 'demographic_info.csv'), encoding='utf_8_sig', index=False)

    # STEP2: 获取各任务完成情况
    task_name = ['02_DDK', '10_PD', '13_DDKWM', '14_PDSS', 'PSQI', 'BDI', 'GAD-7', 'MMSE', 'DST', 'IQCODE']
    pd_status = pd.DataFrame()
    for i_task in task_name:
        tn = i_task.split('_')[-1]
        col = data_dem_speech.columns[data_dem_speech.columns.str.startswith(i_task)]
        valid_id = data_dem_speech.dropna(how='all', subset=col)['id'].to_list()
        i_num = len(valid_id)
        i_num_m = data_dem_speech[(data_dem_speech['sex'] == 1) & (data_dem_speech['id'].isin(valid_id))].shape[0]
        i_num_f = data_dem_speech[(data_dem_speech['sex'] == 0) & (data_dem_speech['id'].isin(valid_id))].shape[0]
        i_age = data_dem_speech[data_dem_speech['id'].isin(valid_id)]['age']
        i_edu = data_dem_speech[data_dem_speech['id'].isin(valid_id)]['edu']
        pd_task = pd.DataFrame({'任务名': [tn], '完成人数 (男/女)': [f"{i_num} ({i_num_m} / {i_num_f})"], 
                                '年龄 (年)': [f"{i_age.mean():.2f} ± {i_age.std():.2f} ({i_age.min()} ~ {i_age.max()})"], 
                                '教育水平 (年)': [f"{i_edu.mean():.2f} ± {i_edu.std():.2f} "
                                             f"({i_edu.min():.2f} ~ {i_edu.max():.2f})"]})
        pd_status = pd.concat([pd_status, pd_task]).reset_index(drop=True)
    pd_status.to_csv(os.path.join(save_dir, 'task_complete_status.csv'), encoding='utf_8_sig', index=False)
    return pd_dem, pd_status


def corr_compare(cog_feat_data: pd.DataFrame, save_dir: Union[str, os.PathLike], wm_thr: int = 5,
                 method: str = 'pearson', padjust: str = 'none') -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    单/双言语特征与认知状况(MMSE/DST/IQCODE)间的相关性比较：对每种任务的言语特征与认知得分的相关系数列表进行配对样本T检验
    :param cog_feat_data: 认知评分与言语特征数据
    :param save_dir: 保存结果路径
    :param wm_thr: 工作记忆次任务的正确数阈值
                   即最少wm_thr个回答正确的被试才被作为待检验的条目（总共10个试次，不少于wm_thr个准确回答才作为有效的双任务范式数据）
    :param method: 计算相关性的方法
    :param padjust: p值校正方法
    :return: 全部认知状态的相关性r值和p值（若校正则为校正后的）; 每种任务的言语特征与认知得分的相关系数相关系数绝对值列表的配对样本T检验结果
    """
    valid_id = cog_feat_data[(cog_feat_data["13_DDKWM-T1-WM-Correct Count"] >= wm_thr) &
                             (cog_feat_data["13_DDKWM-T2-WM-Correct Count"] >= wm_thr)]["id"]
    _cog_feat_data = cog_feat_data.loc[cog_feat_data["id"].isin(valid_id), :]
    pval_stars = {0.001: '***', 0.01: '**', 0.05: '*'}

    def replace_pval(x):
        for key, value in pval_stars.items():
            if x < key:
                return value
        return ''

    cognition = ['MMSE', 'DST', 'IQCODE']
    task_str = {"Single-Task": r"^02_DDK-|^10_PD-", "Dual-Task": r"^13_DDKWM-T\d+-DDK-|^14_PDSS-PD-"}
    df_r_grp, df_p_grp = pd.DataFrame(), pd.DataFrame()
    for t_key, t_val in task_str.items():
        cog_feat = _cog_feat_data[cognition + _cog_feat_data.columns[_cog_feat_data.columns.str.match(t_val)].to_list()]
        if "^13_DDKWM-" in t_val:
            grouped_cols = cog_feat.columns.str.extract(r'^13_DDKWM-T(\d+)-(.*)', expand=True)
            df_mean = cog_feat.groupby(grouped_cols[1].to_list(), axis=1, sort=False).mean()
            df_mean = df_mean.add_prefix("13_DDKWM-")
            cog_feat = cog_feat.drop(columns=cog_feat.columns[cog_feat.columns.str.match("^13_DDKWM-")])
            cog_feat = pd.concat([df_mean, cog_feat], axis=1)
        cog_feat.columns = cog_feat.columns.str.replace(r'^\d{2}_', '', regex=True)\
            .str.replace(r' \(.*\)$', '', regex=True).str.replace(r'^PD-', 'PicD-', regex=True)\
            .str.replace('pause', 'Pause').str.replace('rate', 'Rate').str.replace('duration', 'Duration')\
            .str.replace('regularity', 'Regularity').str.replace(' Mean', '').str.replace(' Num', ' Number')
        _ft = cog_feat.drop(columns=cognition).fillna(cog_feat.drop(columns=cognition).mean())
        cog_feat = pd.concat([cog_feat[cognition], _ft], axis=1)
        acou_feat_rename = cog_feat.drop(columns=cognition).columns.to_list()
        df_r_p_val = cog_feat.rcorr(method=method, upper='pval', decimals=6,
                                    padjust=padjust, stars=False)  # type:pd.DataFrame
        _save_dir = os.path.join(save_dir, t_key)
        if not os.path.exists(_save_dir):
            os.makedirs(_save_dir)
        df_r_p_val.to_csv(os.path.join(_save_dir, "r_p_value.csv"), encoding="utf-8-sig")
        df_pcorr_r, df_pcorr_p = pd.DataFrame(), pd.DataFrame()
        for cog in cognition:
            r_l, p_l, aco_l = [], [], []
            for aco in acou_feat_rename:
                r_l.append(float(df_r_p_val.loc[aco, cog]))
                p_l.append(float(df_r_p_val.loc[cog, aco]))
                aco_l.append(aco)
            df_pcorr_r = pd.concat([df_pcorr_r, pd.DataFrame(data=r_l, index=aco_l, columns=[cog])], axis=1)
            df_pcorr_p = pd.concat([df_pcorr_p, pd.DataFrame(data=p_l, index=aco_l, columns=[cog])], axis=1)
        p_matrix_star = pd.DataFrame(df_pcorr_p.to_numpy()).applymap(replace_pval).to_numpy()
        fig, ax = plt.subplots(figsize=(9, 4), tight_layout=True)
        cbar_kws = {'aspect': 10, "format": "%.2f", "pad": 0.01}
        ax_lo = sns.heatmap(df_pcorr_r.T, annot=p_matrix_star.T, ax=ax, fmt='',
                            annot_kws={'size': 13, 'weight': 'bold', 'color': 'k'}, cmap="coolwarm", cbar_kws=cbar_kws)
        cbar = ax_lo.collections[-1].colorbar
        cbar.ax.tick_params(labelsize=12, length=3)
        cbar.outline.set_visible(False)
        plt.tick_params(bottom=False, left=False)
        ax_lo.set_yticklabels(ax_lo.get_yticklabels(), fontsize=16, fontproperties=font_family, rotation=0)
        ax_lo.set_xticklabels(ax_lo.get_xticklabels(), fontsize=14, fontproperties=font_family, rotation=50,
                              ha="right", rotation_mode="anchor")
        fig_file = os.path.join(_save_dir, 'corr_heatmap.png')
        plt.savefig(fig_file, dpi=600, bbox_inches='tight', pad_inches=0.02)
        plt.savefig(fig_file.replace('.png', '.svg'), format='svg', bbox_inches='tight', pad_inches=0.02)
        plt.savefig(fig_file.replace('.png', '.tif'), dpi=600, bbox_inches='tight', pad_inches=0.02,
                    pil_kwargs={"compression": "tiff_lzw"})
        plt.show()
        plt.close()
        df_r = df_pcorr_r.abs().reset_index().rename(columns={'index': 'features'})
        df_r['group'] = t_key
        df_r_grp = pd.concat([df_r_grp, df_r])
        df_p = df_pcorr_p.abs().reset_index().rename(columns={'index': 'features'})
        df_p['group'] = t_key
        df_p_grp = pd.concat([df_p_grp, df_p])
    df_r_grp.reset_index(inplace=True)
    df_p_grp.reset_index(inplace=True)
    df_r_grp.to_csv(os.path.join(save_dir, "r_abs_all.csv"), encoding="utf-8-sig", index=False)
    df_p_grp.to_csv(os.path.join(save_dir, "p_all.csv"), encoding="utf-8-sig", index=False)
    r_grp = df_r_grp.copy()
    for cog in cognition:
        data_grp = df_r_grp.groupby(['group'])[cog]
        thr_high = data_grp.quantile(0.75) + 1.5 * (data_grp.quantile(0.75) - data_grp.quantile(0.25))
        thr_low = data_grp.quantile(0.25) - 1.5 * (data_grp.quantile(0.75) - data_grp.quantile(0.25))
        for i_grp in thr_high.index:  # 分组异常值赋值为nan并均值填充
            r_grp[cog] = r_grp[cog].mask((r_grp['group'] == i_grp) & ((r_grp[cog] < thr_low[i_grp]) |
                                                                      (r_grp[cog] > thr_high[i_grp])))
    r_grp.fillna(r_grp.groupby(['group']).transform('mean'), inplace=True)
    within, subject = 'group', 'index'
    res = pd.DataFrame()
    for dv in cognition:  # 配对T检验
        for task in ["Dual-Task"]:
            res_t = pg.ttest(r_grp[r_grp[within] == task][dv], r_grp[r_grp[within] == "Single-Task"][dv],
                             paired=True, tail='greater')
            res_t['cognition'] = dv
            res_t['paired T-test'] = f"{task}-Single-Task"
            res = pd.concat([res, res_t]).reset_index(drop=True)
    res = pd.concat([res.iloc[:, -2:], res.iloc[:, :-2]], axis=1)
    res.to_csv(os.path.join(save_dir, "res_ttest.csv"), encoding="utf-8-sig", index=False)
    plt_data = pd.DataFrame({r'$\vert r \vert$': r_grp['MMSE'].to_list() + r_grp['DST'].to_list() + r_grp['IQCODE'].to_list(),
                             'cog': ['MMSE']*len(r_grp['MMSE']) + ['DST']*len(r_grp['DST']) +
                                    ['IQCODE']*len(r_grp['IQCODE']), 'grp': r_grp['group'].to_list() * 3})
    plt.figure(figsize=(6, 6), clear=True, tight_layout=True)
    plot_args = {'data': plt_data, 'x': 'cog', 'y': r'$\vert r \vert$', 'hue': 'grp'}
    ax = sns.boxplot(width=.7, zorder=10, palette="pastel", showfliers=True, capprops={'label': 'cap'}, **plot_args)
    ax = pt.stripplot(width=.7, size=3, dodge=True, ax=ax, zorder=5, **plot_args)
    adjust_box_widths(ax, .7)
    ax.set_ylim(ax.get_ylim()[0], 1.25 * ax.get_ylim()[-1])
    pairs = [(("MMSE", "Single-Task"), ("MMSE", "Dual-Task")), (("DST", "Single-Task"), ("DST", "Dual-Task")),
             (("IQCODE", "Single-Task"), ("IQCODE", "Dual-Task"))]
    annotator = Annotator(ax, pairs, **plot_args)
    plt.xlabel('')
    # plt.ylabel(r'$\vert r \vert$')
    ax.yaxis.label.set_size(18)
    annotator.configure(text_format='star', loc='inside', fontsize=16)
    annotator.set_pvalues(res['p-val']).annotate()
    _handles, _labels = ax.get_legend_handles_labels()
    handles, labels = [], []
    for i_it in range(len(_labels)):
        if not isinstance(_handles[i_it], Line2D):
            handles.append(_handles[i_it])
        if _labels[i_it] != 'cap':
            labels.append(_labels[i_it])
    plt.legend(handles=handles[:len(set(labels))], labels=labels[:len(set(labels))], loc='upper left',
               prop={'family': font_family, 'size': 16}, title=None)
    for sp in plt.gca().spines:
        plt.gca().spines[sp].set_color('k')
        plt.gca().spines[sp].set_linewidth(1)
    plt.xticks(fontproperties=font_family)
    plt.yticks(fontproperties=font_family)
    yticks = ax.yaxis.get_major_ticks()
    yticks[-2].label1.set_visible(False)
    yticks[-2].set_visible(False)
    plt.gca().tick_params(labelsize=14, direction='in', color='k', length=5, width=1)
    plt.tick_params('x', labelsize=18)
    fig_file = os.path.join(save_dir, "ttest.png")
    plt.savefig(fig_file, dpi=600, bbox_inches='tight', pad_inches=0.02)
    plt.savefig(fig_file.replace('.png', '.svg'), format='svg', bbox_inches='tight', pad_inches=0.02)
    plt.savefig(fig_file.replace('.png', '.tif'), dpi=600, bbox_inches='tight', pad_inches=0.02,
                pil_kwargs={"compression": "tiff_lzw"})
    plt.show()
    plt.close()
    return df_r_grp, df_p_grp, res


def reg_compare(cog_feat_data: pd.DataFrame, corr_pval: pd.DataFrame, save_dir: Union[str, os.PathLike],
                wm_thr: int = 5):
    """
    单/双言语特征与认知状况(MMSE/DST/IQCODE)间的回归比较
    :param cog_feat_data: 认知评分与言语特征数据
    :param corr_pval: 单/双言语特征与认知得分间的相关性结果
    :param save_dir: 保存结果路径
    :param wm_thr: 工作记忆次任务的正确数阈值
                   即最少wm_thr个回答正确的被试才被作为待检验的条目（总共10个试次，不少于wm_thr个准确回答才作为有效的双任务范式数据）
    :return: None
    """
    valid_id = cog_feat_data[(cog_feat_data["13_DDKWM-T1-WM-Correct Count"] >= wm_thr) &
                             (cog_feat_data["13_DDKWM-T2-WM-Correct Count"] >= wm_thr)]["id"]
    _cog_feat_data = cog_feat_data.loc[cog_feat_data["id"].isin(valid_id), :]
    cognition = ['MMSE', 'DST', 'IQCODE']
    task_str = {"Single-Task": r"^02_DDK-|^10_PD-", "Dual-Task": r"^13_DDKWM-T\d+-DDK-|^14_PDSS-PD-"}
    plt_data = {'grp': [], 'cog': [], 'aic': [], 'r2': []}
    for cog in cognition:
        cog_data = _cog_feat_data.dropna(subset=[cog])[cog]
        for t_key, t_val in task_str.items():
            speech_data = _cog_feat_data.dropna(subset=[cog]).iloc[:, _cog_feat_data.columns.str.match(t_val)]
            if "^13_DDKWM-" in t_val:
                grouped_cols = speech_data.columns.str.extract(r'^13_DDKWM-T(\d+)-(.*)', expand=True)
                df_mean = speech_data.groupby(grouped_cols[1].to_list(), axis=1, sort=False).mean()
                df_mean = df_mean.add_prefix("13_DDKWM-")
                speech_data = speech_data.drop(columns=speech_data.columns[speech_data.columns.str.match("^13_DDKWM-")])
                speech_data = pd.concat([df_mean, speech_data], axis=1)
            speech_data.columns = speech_data.columns.str.replace(r'^\d{2}_', '', regex=True)\
                .str.replace(r' \(.*\)$', '', regex=True).str.replace(r'^PD-', 'PicD-', regex=True)\
                .str.replace('pause', 'Pause').str.replace('rate', 'Rate').str.replace('duration', 'Duration')\
                .str.replace('regularity', 'Regularity').str.replace(' Mean', '').str.replace(' Num', ' Number')
            thr_high = speech_data.quantile(0.75) + 1.5 * (speech_data.quantile(0.75) - speech_data.quantile(0.25))
            thr_low = speech_data.quantile(0.25) - 1.5 * (speech_data.quantile(0.75) - speech_data.quantile(0.25))
            speech_data = speech_data[(speech_data >= thr_low) & (speech_data <= thr_high)]
            speech_data.fillna(speech_data.mean(), inplace=True)
            valid_col = corr_pval.loc[(corr_pval['group'] == t_key) &
                                      ((corr_pval['MMSE'] <= 0.05) | (corr_pval['DST'] <= 0.05) |
                                       (corr_pval['IQCODE'] <= 0.05)), 'features']
            x_data = speech_data.loc[:, valid_col]
            data_x, data_y = x_data.to_numpy(), cog_data.to_numpy()
            data_x = StandardScaler().fit_transform(data_x)
            reg = SVR(C=3)
            reg.fit(data_x, data_y)
            r2 = reg.score(data_x, data_y)
            y_pred = reg.predict(data_x)
            n = data_x.shape[0]  # 样本数
            sigma2 = mean_squared_error(data_y, y_pred)  # 残差平方和均值，即均方误差
            k = len(reg.support_) + 1  # 模型参数量：支撑向量数量 + 核参数自由度
            aic = n * np.log(sigma2) + 2 * k
            plt_data['grp'].append(t_key)
            plt_data['cog'].append(cog)
            plt_data['aic'].append(aic)
            plt_data['r2'].append(r2)
            # loo = LeaveOneOut()
            # rss = 0
            # for train_index, test_index in loo.split(data_x):
            #     x_train, x_test = data_x[train_index], data_x[test_index]
            #     y_train, y_test = data_y[train_index], data_y[test_index]
            #     scaler = StandardScaler()
            #     x_train = scaler.fit_transform(x_train)
            #     x_test = scaler.transform(x_test)
            #     model = SVR()
            #     model.fit(x_train, y_train)
            #     y_pred = model.predict(x_test)
            #     rss += (y_test[0] - y_pred[0]) ** 2  # 残差平方和
            # n = data_x.shape[0]
            # sigma2 = rss / n  # 样本方差
            # k = len(model.support_) + 1  # 模型参数量：支撑向量数量 + 核参数自由度
            # aic = n * np.log(sigma2) + 2 * k
            # print(f"AIC {t_key}/{cog}: {aic}")
    plt_data = pd.DataFrame(plt_data)
    print(plt_data)
    plt.figure(figsize=(6, 6), tight_layout=True)
    plot_args = {'data': plt_data, 'x': 'cog', 'y': 'aic', 'hue': 'grp'}
    ax_aic = sns.barplot(palette="pastel", **plot_args)
    for i_t in ax_aic.containers:
        for j_cog in i_t:
            ax_aic.text(j_cog.get_x() + j_cog.get_width() / 2, j_cog.get_height(), f'{j_cog.get_height():.2f}',
                        ha='center', va='bottom' if j_cog.get_height() > 0 else 'top', c='k',
                        fontdict={'family': font_family, 'size': 14})
    ax_aic.legend(loc='lower left', prop={'family': font_family, 'size': 16}, title=None)
    ax_aic.set_xlabel('')
    ax_aic.set_ylabel(f'AIC', fontdict={'family': font_family, 'size': 18})
    ax_aic.tick_params(direction='in', color='k', length=5, width=1)
    plt.xticks(fontproperties=font_family, fontsize=18)
    plt.yticks(fontproperties=font_family, fontsize=14)
    ax_r2 = ax_aic.twinx()
    ax_r2.set_ylim(0.0, 1.0)
    for i_cog in range(len(ax_aic.containers[0])):
        x_range, y_val = [], []
        for j_t in range(len(ax_aic.containers)):
            x_range.append(ax_aic.containers[j_t][i_cog].get_x() + ax_aic.containers[j_t][i_cog].get_width() / 2)
            y_val.append(plt_data[plt_data['grp'] == ax_aic.containers[j_t].get_label()]['r2'].tolist()[i_cog])
            if j_t % 2:
                ax_r2.text(x_range[-1], y_val[-1] + 0.02, f'{y_val[-1]:.4f}', ha='center', va='bottom', c='b',
                           fontdict={'family': font_family, 'size': 14})
            else:
                ax_r2.text(x_range[-1], y_val[-1] - 0.02, f'{y_val[-1]:.4f}', ha='center', va='top', c='b',
                           fontdict={'family': font_family, 'size': 14})
        ax_r2.plot(x_range, y_val, 'bo-', lw=2, ms=10)
    ax_r2.set_xlabel('')
    ax_r2.set_ylabel(f'R2', color='b', fontdict={'family': font_family, 'size': 18})
    ax_r2.set_yticks(ax_r2.get_yticks())
    ax_r2.set_yticklabels(ax_r2.get_yticks(), color='b', fontdict={'family': font_family, 'size': 14})
    ax_r2.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
    ax_r2.tick_params(direction='in', colors='b', length=5, width=1)
    for sp in plt.gca().spines:
        plt.gca().spines[sp].set_color('k')
        plt.gca().spines[sp].set_linewidth(1)
    fig_file = os.path.join(save_dir, f"aic.png")
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    plt.savefig(fig_file, dpi=600, bbox_inches='tight', pad_inches=0.02)
    plt.savefig(fig_file.replace('.png', '.svg'), format='svg', bbox_inches='tight', pad_inches=0.02)
    plt.savefig(fig_file.replace('.png', '.tif'), dpi=600, bbox_inches='tight', pad_inches=0.02,
                pil_kwargs={"compression": "tiff_lzw"})
    plt.show()
    plt.close('all')


def age_compare(cog_feat_data: pd.DataFrame, save_dir: Union[str, os.PathLike],
                wm_thr: int = 5, parametric: bool = True):
    """
    绘制每个言语特征在年轻/年老组间的分布与差异
    :param cog_feat_data: 认知评分与言语特征数据
    :param save_dir: 保存结果路径
    :param wm_thr: 工作记忆次任务的正确数阈值
                   即最少wm_thr个回答正确的被试才被作为待检验的条目（总共10个试次，不少于wm_thr个准确回答才作为有效的双任务范式数据）
    :param parametric: 是否使用非参数检验, False时使用非参数检验
    :return: None
    """
    valid_id = cog_feat_data[(cog_feat_data["13_DDKWM-T1-WM-Correct Count"] >= wm_thr) &
                             (cog_feat_data["13_DDKWM-T2-WM-Correct Count"] >= wm_thr)]["id"]
    age_feat_data = cog_feat_data.loc[cog_feat_data["id"].isin(valid_id), :]
    task_str = {"Single-Task": r"^age|^02_DDK-|^10_PD-", "Dual-Task": r"^age|^13_DDKWM-T\d+-DDK-|^14_PDSS-PD-"}
    for t_key, t_val in task_str.items():
        speech_data = age_feat_data.iloc[:, age_feat_data.columns.str.match(t_val)]
        if "^13_DDKWM-" in t_val:
            grouped_cols = speech_data.columns.str.extract(r'^13_DDKWM-T(\d+)-(.*)', expand=True)
            df_mean = speech_data.groupby(grouped_cols[1].to_list(), axis=1, sort=False).mean()
            df_mean = df_mean.add_prefix("13_DDKWM-")
            speech_data = speech_data.drop(columns=speech_data.columns[speech_data.columns.str.match("^13_DDKWM-")])
            speech_data = pd.concat([df_mean, speech_data], axis=1)
        speech_data.columns = speech_data.columns.str.replace(r'^\d{2}_', '', regex=True) \
            .str.replace(r' \(.*\)$', '', regex=True).str.replace(r'^PD-', 'PicD-', regex=True) \
            .str.replace('pause', 'Pause').str.replace('rate', 'Rate').str.replace('duration', 'Duration') \
            .str.replace('regularity', 'Regularity').str.replace(' Mean', '').str.replace(' Num', ' Number')
        thr_high = speech_data.quantile(0.75) + 1.5 * (speech_data.quantile(0.75) - speech_data.quantile(0.25))
        thr_low = speech_data.quantile(0.25) - 1.5 * (speech_data.quantile(0.75) - speech_data.quantile(0.25))
        speech_data = speech_data[(speech_data >= thr_low) & (speech_data <= thr_high)]
        speech_data.fillna(speech_data.mean(), inplace=True)
        age_diff = speech_data.copy()
        col_scale = speech_data.columns.difference(['age'])
        age_diff[col_scale] = StandardScaler().fit_transform(age_diff[col_scale])
        feat_rename = age_diff.drop(columns=['age']).columns
        age_diff['age_group'] = np.where(age_diff['age'] >= 50, 'Older', 'Younger')
        age_diff.drop(columns=['age'], inplace=True)
        age_diff = age_diff.melt(id_vars=['age_group'], var_name='factor_name', value_name='factor_val')
        res = pd.DataFrame()
        for i_ft in feat_rename:
            if parametric:
                res_t = pg.ttest(age_diff[(age_diff['factor_name'] == i_ft) &
                                          (age_diff['age_group'] == 'Younger')]['factor_val'],
                                 age_diff[(age_diff['factor_name'] == i_ft) &
                                          (age_diff['age_group'] == 'Older')]['factor_val'], tail='one-sided')
            else:
                res_t = pg.mwu(age_diff[(age_diff['factor_name'] == i_ft) &
                                        (age_diff['age_group'] == 'Younger')]['factor_val'],
                               age_diff[(age_diff['factor_name'] == i_ft) &
                                        (age_diff['age_group'] == 'Older')]['factor_val'], tail='one-sided')
            res_t['factor_name'] = i_ft
            res_t['T-test'] = 'Younger-Older'
            res = pd.concat([res, res_t]).reset_index(drop=True)
        res = pd.concat([res.iloc[:, -2:], res.iloc[:, :-2]], axis=1)
        plt.figure(figsize=(10, 6), tight_layout=True)
        plot_args = {'data': age_diff, 'x': 'factor_name', 'y': 'factor_val', 'hue': 'age_group'}
        sns.stripplot(dodge=False, alpha=.3, zorder=0, **plot_args)
        ax = sns.pointplot(ci=None, join=False, markers="_", scale=3, **plot_args)
        _handles, _labels = ax.get_legend_handles_labels()
        handles, labels = [], []
        for i_it in range(len(_labels)):
            if not isinstance(_handles[i_it], Line2D):
                handles.append(_handles[i_it])
            if _labels[i_it] != 'cap':
                labels.append(_labels[i_it])
        plt.legend(handles=handles[:len(set(labels))], labels=labels[:len(set(labels))], loc='lower left',
                   prop={'family': font_family, 'size': 16}, frameon=False, title=None,
                   labelspacing=.3, handletextpad=.1)
        pairs = [((i, "Younger"), (i, "Older")) for i in feat_rename]
        annotator = Annotator(ax, pairs, **plot_args)
        annotator.configure(text_format='star', loc='inside', fontsize=15)
        annotator.set_pvalues(res['p-val']).annotate()
        plt.ylim(-4, 4)
        plt.xlabel('')
        plt.ylabel('Standardized Feature Values', fontdict={'family': font_family, 'size': 16})
        plt.gca().tick_params(direction='in', color='k', length=5, width=1)
        plt.xticks(fontproperties=font_family, fontsize=16, rotation=40, ha="right", rotation_mode="anchor")
        plt.yticks(fontproperties=font_family, fontsize=14)
        for sp in plt.gca().spines:
            plt.gca().spines[sp].set_color('k')
            plt.gca().spines[sp].set_linewidth(1)
        fig_file = os.path.join(save_dir, f"{t_key}/age_diff.png")
        if not os.path.exists(os.path.dirname(fig_file)):
            os.makedirs(os.path.dirname(fig_file))
        plt.savefig(fig_file, dpi=600, bbox_inches='tight', pad_inches=0.02)
        plt.savefig(fig_file.replace('.png', '.svg'), format='svg', bbox_inches='tight', pad_inches=0.02)
        plt.savefig(fig_file.replace('.png', '.tif'), dpi=600, bbox_inches='tight', pad_inches=0.02,
                    pil_kwargs={"compression": "tiff_lzw"})
        plt.show()
        plt.close('all')
        res.to_csv(os.path.join(save_dir, f"{t_key}/res_ttest.csv"), encoding="utf-8-sig", index=False)


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


if __name__ == "__main__":
    start_time = datetime.datetime.now()
    print(f"---------- Start Time ({os.path.basename(__file__)}): {start_time.strftime('%Y-%m-%d %H:%M:%S')} ----------")
    current_path = os.path.dirname(os.path.realpath(__file__))
    data_path = os.path.join(current_path, r"data")
    res_path = os.path.join(current_path, r"results")
    # res_path = os.path.join(current_path, r"results/matched")
    data_f = os.path.join(data_path, 'data.csv')

    pd_data = pd.read_csv(data_f)  # STEP1: 获取所有待使用的数据
    # less_than_30 = pd_data[pd_data['age'] < 30].sample(n=30, random_state=rs)
    # pd_data = pd.concat([pd_data[pd_data['age'] >= 30], less_than_30]).sort_index(ignore_index=True)
    # print(pd_data.shape, pd_data[pd_data_match['age'] > 50].shape)  # 98, 49
    run_get_demog_info = False  # STEP2: 获取所有样本的人口统计学信息
    if run_get_demog_info:
        get_demographic_info_task_status(pd_data, min_age=18, age_interval=15, save_dir=res_path,
                                         parametric=True, equ_var=True, padjust=None)
    run_corr_comp = False  # STEP3: 单/双言语特征与认知得分间的相关性比较
    if run_corr_comp:
        r_all, p_all, _ = corr_compare(pd_data, save_dir=os.path.join(res_path, 'corr_compare'),
                                       wm_thr=5, method='spearman')
    else:
        r_all = pd.read_csv(os.path.join(res_path, 'corr_compare/r_abs_all.csv'))
        p_all = pd.read_csv(os.path.join(res_path, 'corr_compare/p_all.csv'))
    run_reg_comp = False  # STEP4: 单/双言语特征回归认知得分比较
    if run_reg_comp:
        reg_compare(pd_data, p_all, save_dir=os.path.join(res_path, 'reg_compare'), wm_thr=5)
    run_age_comp = False  # STEP5: 单/双言语特征在年轻/年老组间比较
    if run_age_comp:
        age_compare(pd_data, save_dir=os.path.join(res_path, 'age_compare'), wm_thr=5, parametric=True)

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
