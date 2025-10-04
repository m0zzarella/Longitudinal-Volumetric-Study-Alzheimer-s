# Longitudinal Volumetric Study for the Progression of Alzheimer’s Disease from Structural MRI

## **2024 IEEE International Conference on Computer Vision and Machine Intelligence (Best Paper Award for Computer Vision)**

Code repository of aforementioned paper for longitudinal processing of tissue volumes in Alzheimer's patients. 

[![paper](https://img.shields.io/badge/Conference-Paper-blue)](https://ieeexplore.ieee.org/document/10782874)
[![arXiv](https://img.shields.io/badge/arXiv-Paper-brightgreen)](https://arxiv.org/abs/2310.05558)
[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)
[![ADNI](https://img.shields.io/badge/ADNI-Database-orange)](https://adni.loni.usc.edu/)

## Abstract

Alzheimer’s Disease (AD) is an irreversible neurode-generative disorder affecting millions of individuals today. The prognosis of the disease solely depends on treating symptoms as they arise and proper caregiving, as there are no current medical preventative treatments apart from newly developing drugs, which can, at most, slow the progression. Thus, early detection of the disease at its most premature state is of paramount importance. We performed a longitudinal study of structural MRI for certain test subjects with AD and temporal data selected randomly from the Alzheimer’s Disease Neuroimaging Initiative (ADNI) database. We implement a robust pipeline to study the data, including modern pre-processing techniques such as spatial image registration, skull stripping, inhomogeneity correction, and tissue segmentation using an unsupervised learning approach (hidden Markov field model) using intensity histogram information. The temporal data across multiple visits is used to study the structural change in volumes of these tissue classes, namely, cerebrospinal fluid (CSF), grey matter (GM), and white matter (WM) as the patients progressed further into the disease. We also analyze the changes in features extracted, with a modified Mann-Kendall statistic to detect changes in volume trends for each patient. A monotonic decrease (85.31\%, 95.54\%) was observed for GM volumes consistent with clinical findings, whereas the CSF counterpart increased monotonically with varying confidence scores (69.85\%, 85.31\%, 95.54\%). Even though it is observed that WM volumes do not follow a majority trend, they exhibit an intolerable shifting from the baseline according to individual prognosis


## Dataset Used
- The Alzheimer’s Disease Neuroimaging Initiative (ADNI) is a longitudinal, multi-center, observational study.
[ADNI](https://adni.loni.usc.edu/)

## Pipeline Functionalities

The repo currently contains the following:

- **DICOM to NIfTI conversion** with dcm2niix(standalone) and dicom2nifti(integrated)
- **N4 bias field correction** using FSL
- **Intensity normalization** with methods limited to z-score, min-max, histogram matching, Nyul
- **Image registration** to Visit1 baseline using FSL FLIRT/FNIRT
- **Skull stripping** with FSL BET and center-of-gravity initialization to improve iterations
- **Tissue segmentation** using FSL FAST (Cerebrospinal fluid, grey matter, white matter)
- **Volume calculation** with longitudinal tracking and delta calculations
- **Trend analysis** using modified Mann-Kendall test for monotonic changes

## Installation

Install RStudio and the following libraries and oackages

### System
- [dcm2niix](https://github.com/rordenlab/dcm2niix) – DICOM → NIfTI conversion  
- [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki) – tools for bias correction, BET, and FAST segmentation  

`FSL` commands such as `bet` and `fast` should be available in `$PATH`.

### R Packages
- oro.dicom  
- oro.nifti  
- neurobase  
- fslr  
- scales  

Install in R:

```R
install.packages(c("oro.dicom", "oro.nifti", "scales"))
remotes::install_github("muschellij2/neurobase")
remotes::install_github("muschellij2/fslr")

```
## Example Folder Structure

Patient folders must contain visit subfolders:

```
Patient1/
  Visit1/
    scan.nii.gz   
  Visit2/
    scan.nii.gz
  Visit3/
    ...
```

## Usage

### Basic dry runs

1. **Convert all DICOMs to NIfTI** (if downloaded files from database already not in NIfTI):

   ```bash
   ./src/scripts/convert_all.sh data/
   ```

3. **Process a single patient**:

    ```bash
   Rscript src/R/03_longitudinal_pipeline.R /path/to/Patient123
   ```

5. **Batch process multiple patients** (uses parallel jobs):

   ```bash
   ./src/scripts/batch_process.sh -d data/ -j 4
   ```


### Detailed Usage 

If you only have DICOMs, run the converter first (see below).

- Step 1: Convert DICOMs to NIfTI

  - Run for a single patient using R script:

    ```bash
    Rscript src/R/01_dicom_to_nifti.R /path/to/Patient123
    ```

    or


    ```bash
    ./src/scripts/convert_single_patient.sh data/Patient123
    ```
    
  - Convert all patients simultaneously:

    ```bash
    ./src/scripts/convert_all.sh data/
    ```

    This creates .nii.gz files inside each Visit folder.

- Step 2: Longitudinal Processing

  - Run the full pipeline (includes intensity normalization and registration also):

  ```bash
  Rscript src/R/03_longitudinal_pipeline.R /path/to/Patient123
  ```

- Step 3: Standalone Runs of Individual Files

  - Intensity normalization:

  ```bash
  Rscript src/R/04_intensity_normalization.R input.nii.gz output.nii.gz [method] [mask.nii.gz] [reference.nii.gz]
  ```

  Methods: `zscore`, `minmax`, `histogram_match`, `nyul` (default: `zscore`)

  - Image registration:

  ```bash
  Rscript src/R/05_registration.R moving.nii.gz fixed.nii.gz output.nii.gz [linear|nonlinear|both] [cost_function]
  ```

  Cost functions: `corratio`, `mutualinfo`, `normmi`, `normcorr`, `leastsq` (default: `corratio`)

  - Outputs stored as below:

  ```
  results/Patient123/
    Visit1/
      Visit1_N4.nii.gz                    #bias-corrected
      Visit1_N4_norm.nii.gz               #intensity normalized (if needed)
      Visit1_BET.nii.gz                   #skull stripped (baseline)
      Visit1_BET_pve_{0,1,2}.nii.gz       #CSF/GM/WM segmentations
      Visit1_volumes.csv
      Visit1_qc_overlay.png               #tissue segmentation QC
    Visit2/
      Visit2_N4.nii.gz                    #bias-corrected
      Visit2_N4_norm.nii.gz               #intensity normalized
      Visit2_BET.nii.gz                   #skull stripped
      Visit2_BET_reg.nii.gz               #registered to Visit1
      Visit2_registration_quality.csv     #registration metrics
      Visit2_registration_qc.png          #registration QC overlay
      Visit2_BET_pve_{0,1,2}.nii.gz       #CSF/GM/WM segmentations
      Visit2_volumes.csv
      Visit2_qc_overlay.png               #tissue segmentation QC
    Visit3/
      ...
    Patient123_longitudinal_volumes.csv  #per-visit volumes
    Patient123_longitudinal_deltas.csv   #changes from baseline
  ```


### Batch Processing

To process multiple patients in parallel:

```bash

./src/scripts/batch_process.sh -d data/ -j 8 #4 recommended as this takes some time

#skip conversion of dicom to nifti if already done 
./src/scripts/batch_process.sh -d data/ --skip-conversion

./src/scripts/batch_process.sh -d data/ --skip-pipeline
```


### Makefile 

Use make commands as given below for setup:

```bash
#environment setup
make setup

#single patient
make convert PATIENT=Patient123
make process PATIENT=Patient123
make full PATIENT=Patient123        #convert + process

#batch operations
make batch DATA_DIR=data/
make convert-all DATA_DIR=data/

#quality control and maintenance
make qc
make clean                          #dry run cleanup
make clean-force                    #actual cleanup
make status                         #show pipeline status
make report                         #generate QC report

#statistical analysis
make mann-kendall                   #analyze all patients
make mann-kendall-patient PATIENT=Patient123  #analyze specific patient
make mann-kendall-check             #check dependencies
make mann-kendall-preview           #preview analysis (dry run)
```

## Mann-Kendall Trend Analysis

A statistical trend analysis module is also included which implements a modified Mann-Kendall test to detect monotonic trends in brain tissue volume changes over time. This analysis helps identify whether patients show significant volume increases or decreases across different visits.

### Methodology

The analysis uses a modified Mann-Kendall test with the following highlighted features:

1. **Modified Sign Function**: Uses a 1mL threshold to ignore minor fluctuations:
   - `+1` if volume difference > 1mL (meaningful increase)
   - `-1` if volume difference < -1mL (meaningful decrease)  
   - `0` if difference is within [-1mL, +1mL] (no meaningful change)

2. **Statistical Computation**:
   - Mann-Kendall statistic: `S = Σ sgn(x_j - x_k)` for all pairs j > k
   - Variance: `Var(S) = n(n-1)(2n+5)/18`
   - Z-score: `Z = (S-1)/√Var(S)` if S>0, `Z = (S+1)/√Var(S)` if S<0, `Z = 0` if S=0

3. **Trend Interpretation**:
   - **Positive Z**: Volume tends to increase over visits
   - **Negative Z**: Volume tends to decrease over visits
   - **|Z| ≥ 1.96**: Statistically significant at α = 0.05 (two-tailed)

### Usage

#### Run Analysis on All Patients
```bash
make mann-kendall
# or directly:
./src/scripts/mann_kendall_analysis.sh results/
```

#### Analyze Specific Patient
```bash
make mann-kendall-patient PATIENT=Patient123
# or directly:
./src/scripts/mann_kendall_analysis.sh results/ Patient123
```

#### Check Dependencies and Preview
```bash
make mann-kendall-check    
make mann-kendall-preview  
```

### Output Files

The analysis generates:

1. **`mann_kendall_analysis_results.csv`**: Detailed results for each patient and tissue type
   - Patient_ID, Tissue_Type, N_Visits
   - S_Statistic, Variance_S, Z_Score, P_Value
   - Confidence_Percentage, Kendall_Tau, Trend_Interpretation
   - Visit_Sequence, Volume_Sequence

2. **`mann_kendall_summary.txt`**: Explainable and interpretable summary with:
   - Overall statistics (total analyses, significant trends)
   - Breakdown by tissue type (CSF, GM, WM)
   - Patients with multiple significant trends

### Example Output

```
=== Analyzing patient: Patient123 ===
Found 4 visits: Visit1, Visit2, Visit3, Visit4
CSF: S=6, Z=2.236, p=0.0253, confidence=98.74%
GM:  S=-4, Z=-1.342, p=0.1797, confidence=91.01%
WM:  S=-2, Z=-0.671, p=0.5023, confidence=74.89%
```

### Clinical Interpretation

- **CSF Volume Increases**: This is associated with brain atrophy and ventricular enlargement, fills up spaces left from brain tissue atrophy
- **GM Volume Decreases**: Consistent with neurodegeneration in AD progression  
- **WM Volume Changes**: Variable patterns depending on disease stage and individual factors, consistent intolerable shift from baseline - can be due to legion formation or hyperintensive white matter disease correlated with dementia




Paper: 

```
@inproceedings{10782874,
  author={Sanyal, Prayas and Mukherjee, Srinjay and Das, Arkapravo and Sen, Anindya},
  booktitle={2024 IEEE International Conference on Computer Vision and Machine Intelligence (CVMI)}, 
  title={Longitudinal Volumetric Study for the Progression of Alzheimer's Disease from Structural MRI}, 
  year={2024},
  pages={1-8},
  doi={10.1109/CVMI61877.2024.10782874}}
```

## Acknowledgements


- **Alzheimer's Disease Neuroimaging Initiative (ADNI)** for providing the longitudinal neuroimaging dataset used in this study. Data collection was made possible for the access to the ADNI database after approval.

- **Neurohacking** ([https://github.com/muschellij2/Neurohacking](https://github.com/muschellij2/Neurohacking)) - Extensive code adaptations were made from this comprehensive neuroimaging analysis repository from John Hopkins University available in coursera.

- **fslr** ([https://github.com/muschellij2/fslr](https://github.com/muschellij2/fslr)) - Core FSL interface functions and implementations were adapted from this R package for neuroimaging analysis.
