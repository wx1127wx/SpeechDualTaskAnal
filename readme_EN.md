# SpeechDualTaskAnal
Python code for speech–cognition dual-task analysis (speech, pupil, t-fMRI)

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg?style=flat)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.8%2B-orange.svg)](https://www.python.org/)
[![openSMILE](https://img.shields.io/badge/openSMILE-2.x-brightgreen.svg)](https://www.audeering.com/opensmile/)
[![Librosa](https://img.shields.io/badge/Librosa-0.9%2B-green.svg)](https://github.com/librosa/librosa)

> If this repo helps, a ⭐️ is appreciated.

## Introduction
This repository provides companion code for a speech-based dual-task paradigm (e.g., **DDK-WM / PDSS**) used in our study. It includes:
- **Acoustic features** via Librosa and/or openSMILE: prosodic (duration, energy, ZCR, F0), spectral (MFCCs), and voice-quality (formants, jitter, shimmer).
- **Pupillometry (pupil/gaze/fixation)** aligned to task windows: export, blink handling, interpolation, filtering, and aggregated statistics.
- **t-fMRI utilities**: generation of condition timing files and automated AFNI commands (e.g., `3dttest++`, `3dANOVA*`) for group-level analyses and connectivity workflows.

## Components (overview)
- **Part 1 | Speech × Eye-tracking**: behavioral statistics/visualization, pupil data export & preprocessing, task-window slicing and merging.  
- **Part 2 | t-fMRI & AFNI**: timing file creation from event tables; automatic generation of 3dttest++ / 3dANOVA commands.  
- **Part 3 | Integrated Analyses**: demographics summary, correlation and regression comparisons, between-group differences.

## References & Acknowledgements (code/methods sources)
- **openSMILE & canonical feature sets** (eGeMAPS, ComParE, IS09)  
  - Eyben et al., *IEEE TAFFC*, 2015 (GeMAPS)  
  - Schuller et al., *INTERSPEECH 2009/2013/2016* (IS09/ComParE)
- **Librosa**: McFee et al., *Librosa: Audio and Music Signal Analysis in Python*  
- **Pupil Labs**: export fields/timestamps and synchronization practices  
- **AFNI**: official docs and example scripts for `3dttest++`, `3dANOVA*`, gPPI-style analyses  
- **Statistics & Visualization**: `pingouin` (RM-ANOVA/effect sizes/post-hoc), `statannotations` (significance marks), `ptitprince` (rainclouds), `matplotlib` / `seaborn`

We thank the maintainers and contributors of the above open-source projects. Please follow their LICENSE terms when citing and reusing.

## Dependencies (example)
- Python ≥ 3.8  
- `numpy`, `pandas`, `scipy`, `matplotlib`, `seaborn`  
- `pingouin`, `statannotations`, `ptitprince`  
- Optional: `librosa`, `opensmile` (or openSMILE binaries/configs)

## License
MIT License. When using or modifying this code, please keep the copyright and license notice and acknowledge relevant sources in your publications.
