# SpeechDualTaskAnal
Python code for speech–cognition dual-task analysis (speech, pupil, t-fMRI)  

# README: Full Reproduction Pipeline (Eye → Imaging → Statistics)

---

## Part 1: Eye-tracking Data Processing

### Input Data (CSV)
This step uses the following files as input:

**Blink-related**
- `pupil_blink.csv`: Blink-related metrics  
- `blink_event_detail.csv`: Detailed blink event table  
- `blink_data_use_raw.csv`: Integrated blink table used for later matching and analysis  

**Pupil / HiPA metrics**
- `pupil_hipa.csv`: Pupil HiPA metrics  
- `hipa_all.csv`: HiPA summary table  

**Time-window / task matching**
- `data_match30.csv`: Matching table aligning eye metrics with experimental task time windows  

**Subject / task information**
- `subinfo_speech.csv`: Subject metadata and task condition information  

### Processing Logic (High-level)
1. Read blink + pupil + HiPA CSV files  
2. Align metrics to experimental task time windows using `data_match30.csv`  
3. Merge with `subinfo_speech.csv` to append subject and task-condition metadata  
4. Export the integrated eye-tracking analysis table  

### Output
- **Integrated eye-tracking dataset** containing blink metrics, pupil/HiPA metrics, time-window alignment, and subject/task information  

---

## Part 2: Imaging Processing (Four Steps)

**timing files → subject-level preprocessing/first-level model → group-level t-test → ANOVA/PPI and visualization**

---

### Step 1: Generate Timing Files

**Input**
- `timing_file/`: Raw information/intermediate tables used to generate timing files  

**Code**
- `creat_timing_file.py`

**Run**
- Execute `creat_timing_file.py`

**Output**
- Timing files (event onset/offset and condition information), used in Step 2 first-level modeling  

---

### Step 2: Subject-level Processing (Preprocessing + First-level Model)

**Input**
- Timing files generated in Step 1  
- Raw imaging data for each subject  

**Code**
- `afni_proc.sh`  
- `proc.default_subid`  
- `proc_run.sh`  

**Run**
1. Configure/execute `afni_proc.sh`  
2. Run `proc_run.sh` to execute the full subject-level workflow  
3. `proc.default_subid` is used for default subject IDs / batch processing of multiple subjects  

**Output**
- Subject-level preprocessing results and first-level model outputs (standard AFNI outputs: preprocessed time series, design matrices, regression outputs, beta/t-stat maps, etc.)  
- Subject-level contrast/statistical files for group-level analysis (Step 3/Step 4)  

---

### Step 3: Group-level Statistics (t-test) + Visualization Export

**Input**
- Subject-level outputs from Step 2 (contrast/statistical files)  
- `roi_mask/`: Group-level mask/ROI  

**Code**
- `3dttest_autocmd.py`  
- `group_analysis.ttest`  
- `ttest_run.sh`  
- `cmd.txt`  
- `batch_save_surfvol_figs.sh`  
- `BrainNetOption.mat`  

**Run**
1. Run `3dttest_autocmd.py` to automatically generate/organize t-test commands  
2. Run `ttest_run.sh` to execute the group-level t-test  
3. Run `batch_save_surfvol_figs.sh` to batch export surface/volume visualization figures  

**Output**
- `results/`: Group-level t-test statistical results (t-maps, thresholds, clusters, etc.)  
- Visualization figures  
- `cmd.txt`: Command log  

---

### Step 4: Further Group-level Analysis (ANOVA / PPI) + Visualization Export

**Input**
- Results from Step 2  
- `roi_mask/`  

**Code**
- `3danova_autocmd.py`  
- `s1.ppi.sh`  
- `cmd.ppi.2.make.regs`  
- `proc.3dd.ppi.post.full`  
- `s2.3danova.sh`  
- `cmd.txt`  
- `s3.batch_save_surfvol_figs.sh`  
- `BrainNetOption.mat`  

**Run (in filename order)**
1. Run `s1.ppi.sh` (PPI-related processing: generate regressors, etc.)  
2. Run `s2.3danova.sh` (execute group-level 3dANOVA statistics)  
3. Run `s3.batch_save_surfvol_figs.sh` (batch export visualization figures)  
4. `3danova_autocmd.py` automatically generates/organizes commands  

**Output**
- `results/`: Group-level ANOVA/PPI statistical results  
- Visualization figures (surface/volume)  
- `cmd.txt`: Command log  

---

## Part 3: Statistical Analysis and Visualization

### Goal
Use a master script to perform statistical analysis and generate figures based on the integrated dataset, including:

- Demographic and task-completion summary  
- Single-task vs dual-task: correlation between speech features and cognitive scales (heatmaps + paired tests)  
- Single vs dual-task: regression performance comparison (SVR + AIC)  
- Young vs old group: speech-feature differences and visualization  

### Input
- `data.csv`

### Code
- `main.py`

### Run
- Execute `python main.py`

### Output
The script generates tables and figures used for paper reproduction (output directory and filenames depend on script settings), typically including:

- Demographic and task-completion summary table (CSV)  
- Single/dual-task correlation matrices and heatmaps (png/svg/tif)  
- Statistical test results comparing correlation strength (CSV) and corresponding plots  
- Regression comparison (AIC/model performance) plots  
- Age-group difference plots and statistical test tables  
- `finished.txt`: Records script completion time and runtime (reproducibility log)  

---

## Summary: Recommended Reproduction Order

1. **Part 1 (Eye)** → Generate integrated eye-tracking dataset  
2. **Part 2 (Imaging)** → Generate subject-level + group-level statistical results and figures  
3. **Part 3 (Statistics)** → Read `data.csv` and output final statistical tables and figures (main results of the paper)  
