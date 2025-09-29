# SpeechDualTaskAnal
Python code for speech–cognition dual-task analysis (speech, pupil, t-fMRI)

## 简介
本仓库为论文配套代码，面向 **DDK-WM / PDSS** 等言语-工作记忆双任务范式。涵盖：
- **语音声学特征**：基于 Librosa 与/或 openSMILE，提取韵律（时长、能量、过零率、F0）、谱相关（MFCC）与音质（共振峰、jitter、shimmer）等；
- **眼动数据（瞳孔/凝视/注视）**：对齐任务窗口的导出与预处理（眨眼处理、插值、滤波、聚合统计）；
- **t-fMRI**：由事件表生成条件 timing 文件，并自动生成 AFNI 统计命令（如 `3dttest++`、`3dANOVA*`）以复现组层面分析与连通性流程。

> 如果对你有帮助，欢迎点个 Star :)

## 组成（概览）
- **Part 1｜语音 × 眼动**：行为统计、可视化、瞳孔数据导出与预处理、任务窗口切片与合并  
- **Part 2｜t-fMRI & AFNI**：基于事件表创建 timing 文件，自动生成 3dttest++ / 3dANOVA 系列命令  
- **Part 3｜综合分析**：人口学汇总、相关与回归比较、组间差异等常见统计结果

## 参考与致谢（代码/方法来源）
- **openSMILE 与经典特征集**（eGeMAPS、ComParE、IS09）  
  - Eyben et al., *IEEE TAFFC*, 2015（GeMAPS）  
  - Schuller et al., *INTERSPEECH 2009/2013/2016*（IS09/ComParE）
- **Librosa**：McFee et al., *Librosa: Audio and Music Signal Analysis in Python*  
- **Pupil Labs**：导出字段/时间戳与同步思路参考其文档与社区实践  
- **AFNI**：`3dttest++`、`3dANOVA*`、gPPI 等官方手册与示例脚本  
- **统计与可视化**：`pingouin`（RM-ANOVA/效应量/事后比较）、`statannotations`（显著性标注）、`ptitprince`（raincloud 图）、`matplotlib`/`seaborn`

感谢上述开源项目与作者的贡献。请遵循各自 LICENSE 进行引用与再利用。

## 许可
GPL-3.0。使用或修改本仓库代码时，请保留版权与许可声明，并在论文/报告的方法部分致谢与引用相关来源。
