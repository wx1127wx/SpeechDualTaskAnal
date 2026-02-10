# SpeechDualTaskAnal
Python code for speech–cognition dual-task analysis (speech, pupil, t-fMRI)
ReadMe Language | [English readme](./readme_EN.md) | 

# README：完整复现流程（按顺序：眼动 → 影像 → 统计）

本文档仅说明 **“使用什么数据 → 运行什么代码 → 得到什么结果”**。不约束你的目录结构；读者只需根据脚本内的路径设置，将对应数据放到脚本能读到的位置即可。

---

## Part 1：眼动数据处理（Blink / Pupil / HiPA 数据整合）

### 目标
将眼动相关的多个 CSV 中间表进行匹配与融合，得到后续统计/影像分析可直接使用的眼动指标数据。

### 输入数据（CSV）
本步骤使用以下文件作为输入（文件名与内容含义如下）：

**眨眼相关**
- `pupil_blink.csv`：眨眼相关指标（事件级/汇总级，取决于你的生成方式）
- `blink_event_detail.csv`：眨眼事件详细表（起止/持续等事件信息）
- `blink_data_use_raw.csv`：用于后续匹配与分析的眨眼整合表

**瞳孔 / HiPA 指标**
- `pupil_hipa.csv`：瞳孔 HiPA（或瞳孔唤醒/负荷相关）指标
- `hipa_all.csv`：HiPA 汇总表（跨被试/跨条件整合）

**时间窗/任务匹配**
- `data_match30.csv`：将眼动指标与实验任务时间窗口对齐的匹配表（例如 30 秒窗口）

**被试/任务信息**
- `subinfo_speech.csv`：被试信息与任务条件信息（用于补充元信息与分组变量）

### 使用代码
- （你后续提供该步骤实际脚本名/主入口后，我会把脚本名补到此处并写成“如何运行”一句话版本。）

### 处理逻辑（高层）
1. 读取 blink + pupil + HiPA 相关 CSV
2. 依据 `data_match30.csv` 将指标与任务时间窗口/区段对齐
3. 与 `subinfo_speech.csv` 合并，补充被试与任务条件信息
4. 导出整合后的眼动分析表

### 输出结果
- **眼动整合数据表（analysis-ready）**：包含眨眼指标、瞳孔/HiPA 指标、时间窗匹配信息及被试/任务信息  
  （输出文件名以你的脚本实际导出为准；建议在论文复现时保留该输出表作为后续输入之一。）

---

## Part 2：影像处理（AFNI 工作流，四步）

### 目标
完成 fMRI 的 AFNI 分析流程：  
**timing files → 个体层面预处理/一阶模型 → 组水平 t-test → 进一步 ANOVA/PPI 与可视化**

---

### Step 1：生成 timing files

**输入**
- `timing_file/`：生成时序文件所需的原始信息/中间表（内容由你的实验范式决定）

**代码**
- `creat_timing_file.py`

**运行**
- 运行 `creat_timing_file.py`

**输出**
- timing files（事件起止时间/条件信息等），供 Step 2 个体层面模型调用

---

### Step 2：个体层面处理（预处理 + 一阶模型）

**输入**
- Step 1 生成的 timing files
- 每个被试的原始影像数据（路径在脚本中配置）

**代码**
- `afni_proc.sh`
- `proc.default_subid`
- `proc_run.sh`

**运行（推荐顺序）**
1. 执行/配置 `afni_proc.sh`
2. 运行 `proc_run.sh`（执行完整个体层面流程）
3. `proc.default_subid` 用于默认被试编号/批处理配置（如需要批处理多个被试）

**输出**
- 每个被试的预处理结果与一阶模型结果（AFNI 标准输出：预处理后的时序、设计矩阵/回归输出、beta/tstat 等）
- 供 Step 3/Step 4 组水平分析使用的个体层面对比/统计文件

---

### Step 3：组水平统计（t-test） + 可视化导出

**输入**
- Step 2 产生的个体层面结果（被试层面的对比/统计文件）
- `roi_mask/`（可选）：组水平 mask/ROI

**代码**
- `3dttest_autocmd.py`
- `group_analysis.ttest`
- `ttest_run.sh`
- `cmd.txt`
- `batch_save_surfvol_figs.sh`
- `BrainNetOption.mat`

**运行**
1. 运行 `3dttest_autocmd.py` 自动生成/整理 t-test 命令（常写入 `cmd.txt` 或生成可运行脚本）
2. 运行 `ttest_run.sh` 执行组水平 t-test
3. 运行 `batch_save_surfvol_figs.sh` 批量导出 surface/volume 可视化图（如需）

**输出**
- `results/`：组水平 t-test 统计结果（t-map/阈值/聚类等）
- 可视化图片（surface/volume）
- `cmd.txt`：命令记录（建议保留作为复现证据链）

---

### Step 4：进一步组水平分析（ANOVA / PPI） + 可视化导出

**输入**
- Step 2（必要时也包括 Step 3）的结果文件
- `roi_mask/`（可选）

**代码**
- `3danova_autocmd.py`
- `s1.ppi.sh`
- `cmd.ppi.2.make.regs`
- `proc.3dd.ppi.post.full`
- `s2.3danova.sh`
- `cmd.txt`
- `s3.batch_save_surfvol_figs.sh`
- `BrainNetOption.mat`

**运行（按文件名顺序）**
1. 运行 `s1.ppi.sh`（PPI 相关步骤：生成回归量/寄存器等）
2. 运行 `s2.3danova.sh`（执行 3dANOVA 组水平统计）
3. 运行 `s3.batch_save_surfvol_figs.sh`（批量导出可视化）
4. `3danova_autocmd.py` 用于自动生成/整理命令（常写入 `cmd.txt`）

**输出**
- `results/`：ANOVA/PPI 相关的组水平统计结果
- 可视化图片（surface/volume）
- `cmd.txt`：命令记录（复现证据链）

---

## Part 3：统计分析与作图（Python 一体化脚本）

### 目标
基于一个总控脚本对整合后的数据进行统计与可视化输出，包括：
- 人口统计学与任务完成情况汇总
- 单任务 vs 双任务：语音特征与认知量表相关性对比（热图 + 配对检验）
- 单/双任务：回归能力对比（SVR + AIC）
- 年轻/年老组：语音特征差异分析与作图

### 输入数据
- `data.csv`（本部分唯一必需输入）

### 使用代码
- `main.py`

### 运行方式
- 运行 `python main.py`  
  （脚本内部会读取 `data.csv` 并将所有结果写入脚本指定的 `results/` 目录；路径以脚本内设置为准。）

### 输出结果
脚本会生成用于论文复现的一系列表格与图片（输出目录与文件名以脚本实际写入为准），典型包括：
- 人口统计学与任务完成情况汇总表（CSV）
- 单/双任务相关矩阵与热图（png/svg/tif）
- 单/双任务相关强度比较的统计检验结果（CSV）及图
- 回归比较（AIC/模型表现）图
- 年龄组差异图与对应检验结果表
- `finished.txt`：记录脚本完成时间与耗时（复现实验日志）

---

## 总结：复现顺序

1. **Part 1（眼动）**：得到眼动整合表（供后续融合/统计/建模）
2. **Part 2（影像）**：得到个体层面结果 + 组水平统计结果与图
3. **Part 3（统计）**：读取 `data.csv` 输出统计表与图（论文主结果图表的主要来源）

