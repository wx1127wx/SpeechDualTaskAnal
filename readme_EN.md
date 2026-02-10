# SpeechDualTaskAnal
Python code for speech–cognition dual-task analysis (speech, pupil, t-fMRI)  
ReadMe Language | [English readme](./readme_EN.md) | 

# README: Full Reproducible Pipeline (in order: Eye-tracking → Neuroimaging → Statistics)

This document only describes **“what data to use → what code to run → what results are obtained”**.  
No directory structure is required; readers only need to ensure the scripts can access the corresponding data according to the paths defined inside the scripts.

---

## Part 1: Eye-tracking Data Processing (Blink / Pupil / HiPA Integration)

### Goal
Integrate multiple eye-tracking intermediate CSV tables to generate eye-movement metrics that can be directly used for subsequent statistical and neuroimaging analyses.

### Input Data (CSV)
This step uses the following files as input (filenames and meanings are listed below):

**Blink-related**
- `pupil_blink.csv`: Blink-related metrics (event-level or summary-level depending on your processing method)
- `blink_event_detail.csv`: Detailed blink event table (onset/offset/duration and event information)
- `blink_data_use_raw.csv`: Integrated blink table used for downstream matching and analysis

**Pupil / HiPA metrics**
- `pupil_hipa.csv`: Pupil-derived HiPA (or pupil arousal/workload-related) metrics
- `hipa_all.csv`: Aggregated HiPA table (across subjects/conditions)

**Time window / task alignment**
- `data_match30.csv`: Mapping table aligning eye-tracking metrics with experimental time windows (e.g., 30-second windows)

**Subject / task metadata**
- `subinfo_speech.csv`: Subject and task information (used to supplement metadata and grouping variables)

### Code
- (After you provide the actual script name/entry for this step, it can be inserted here with a one-line “how to run” instruction.)

### Processing Logic (High-level)
1. Load blink + pupil + HiPA related CSV tables
2. Use `data_match30.csv` to align metrics with task time windows/segments
3. Merge with `subinfo_speech.csv` to attach subject and task metadata
4. Export the integrated eye-tracking analysis table

### Output
- **Integrated eye-tracking dataset (analysis-ready)**: contains blink metrics, pupil/HiPA metrics, time-window alignment, and subject/task information  
  (Output filename depends on your script; it is recommended to retain this table as a key downstream input for reproducibility.)

---

## Part 2: Neuroimaging Processing (AFNI Workflow, Four Steps)

### Goal
Complete the AFNI fMRI analysis pipeline:  
**timing files → first-level preprocessing/modeling → group-level t-test → extended ANOVA/PPI analysis with visualization**

---

### Step 1: Generate Timing Files

**Input**
- `timing_file/`: Raw/intermediate data required to generate timing files (depends on your experimental design)

**Code**
- `creat_timing_file.py`

**Run**
- Execute `creat_timing_file.py`

**Output**
- Timing files (event onset/duration/condition information) for first-level modeling in Step 2

---

### Step 2: First-level Processing (Preprocessing + First-level Model)

**Input**
- Timing files generated in Step 1
- Raw fMRI data for each subject (paths defined in scripts)

**Code**
- `afni_proc.sh`
- `proc.default_subid`
- `proc_run.sh`

**Run (recommended order)**
1. Execute/configure `afni_proc.sh`
2. Run `proc_run.sh` (executes the full first-level pipeline)
3. `proc.default_subid` is used for default subject ID or batch configuration (when processing multiple subjects)

**Output**
- Preprocessed time series and first-level GLM results for each subject (AFNI standard outputs: design matrix, regression outputs, beta/tstat maps, etc.)
- Subject-level contrast/statistical maps for use in Step 3 and Step 4

---

### Step 3: Group-level Statistics (t-test) + Visualization

**Input**
- Subject-level results from Step 2
- `roi_mask/` (optional): mask/ROI for group analysis

**Code**
- `3dttest_autocmd.py`
- `group_analysis.ttest`
- `ttest_run.sh`
- `cmd.txt`
- `batch_save_surfvol_figs.sh`
- `BrainNetOption.mat`

**Run**
1. Run `3dttest_autocmd.py` to automatically generate/organize t-test commands (often written to `cmd.txt`)
2. Run `ttest_run.sh` to perform group-level t-test
3. Run `batch_save_surfvol_figs.sh` to export surface/volume visualization (if needed)

**Output**
- `results/`: group-level t-test statistical results (t-maps/thresholded maps/clusters, etc.)
- Visualization figures (surface/volume)
- `cmd.txt`: command log (recommended to retain as a reproducibility record)

---

### Step 4: Extended Group-level Analysis (ANOVA / PPI) + Visualization

**Input**
- Results from Step 2 (and Step 3 if required)
- `roi_mask/` (optional)

**Code**
- `3danova_autocmd.py`
- `s1.ppi.sh`
- `cmd.ppi.2.make.regs`
- `proc.3dd.ppi.post.full`
- `s2.3danova.sh`
- `cmd.txt`
- `s3.batch_save_surfvol_figs.sh`
- `BrainNetOption.mat`

**Run (in order)**
1. Run `s1.ppi.sh` (PPI-related processing: generate regressors, etc.)
2. Run `s2.3danova.sh` (execute group-level 3dANOVA)
3. Run `s3.batch_save_surfvol_figs.sh` (export visualization)
4. `3danova_autocmd.py` is used to automatically generate/organize commands (usually written to `cmd.txt`)

**Output**
- `results/`: ANOVA/PPI group-level statistical results
- Visualization figures (surface/volume)
- `cmd.txt`: command record (reproducibility evidence)

---

## Part 3: Statistical Analysis and Visualization (Python Unified Script)

### Goal
Run a unified statistical pipeline using one main script to produce:
- Demographic and task completion summary
- Single-task vs Dual-task correlation comparison between speech features and cognitive scores (heatmaps + paired tests)
- Regression performance comparison (SVR + AIC)
- Age-group (Young vs Older) speech feature comparison and visualization

### Input Data
- `data.csv` (the only required input for this step)

### Code
- `main.py`

### Run
- Execute `python main.py`  
  (The script reads `data.csv` and writes all results into the `results/` directory defined inside the script.)

### Output
The script generates tables and figures used to reproduce the manuscript results, including:
- Demographic and task completion summary tables (CSV)
- Single/Dual-task correlation matrices and heatmaps (png/svg/tif)
- Statistical comparison of correlation strength (CSV + plots)
- Regression comparison (AIC/model performance) figures
- Age-group comparison figures and statistical result tables
- `finished.txt`: execution time and completion log (reproducibility record)

---

## Summary: Recommended Reproduction Order

1. **Part 1 (Eye-tracking)**: generate integrated eye-tracking table (for downstream fusion/statistics/modeling)
2. **Part 2 (Neuroimaging)**: generate subject-level and group-level statistical results and figures
3. **Part 3 (Statistics)**: run `data.csv` to generate statistical tables and figures (primary manuscript outputs)

---

## References and Acknowledgements (Code / Methods)

- **openSMILE and classic feature sets** (eGeMAPS, ComParE, IS09)  
  - Eyben et al., *IEEE TAFFC*, 2015 (GeMAPS)  
  - Schuller et al., *INTERSPEECH 2009/2013/2016* (IS09/ComParE)
- **Librosa**: McFee et al., *Librosa: Audio and Music Signal Analysis in Python*  
- **Pupil Labs**: field definitions, timestamps, and synchronization concepts follow official documentation and community practices  
- **AFNI**: `3dttest++`, `3dANOVA*`, gPPI and related official manuals and example scripts  
- **Statistics and visualization**: `pingouin` (RM-ANOVA/effect size/post-hoc), `statannotations` (significance annotation), `ptitprince` (raincloud plots), `matplotlib`/`seaborn`

We thank the developers of these open-source projects. Please follow their respective LICENSE terms when using or redistributing.

---

## Dependencies (Example)

- Python ≥ 3.8  
- `numpy`, `pandas`, `scipy`, `matplotlib`, `seaborn`  
- `pingouin`, `statannotations`, `ptitprince`  
- Optional: `librosa`, `opensmile` (or openSMILE executable/config)

---

## License

MIT License. When using or modifying this repository, please retain copyright and license notices, and acknowledge relevant sources in the Methods section of papers/reports.
