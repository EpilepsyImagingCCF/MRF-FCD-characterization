# MRF-FCD-characterization
Code for publication on Ann Neurol: Multiparametric Characterization of Focal Cortical Dysplasia Using 3D MR Fingerprinting

This repository accompanies the research paper:

Su, T.-Y., Choi, J.Y., Hu, S., Wang, X., Blümcke, I., Chiprean, K., Krishnan, B., Ding, Z., Sakaie, K., Murakami, H., Alexopoulos, A.V., Najm, I., Jones, S.E., Ma, D. and Wang, Z.I.  
*Multiparametric Characterization of Focal Cortical Dysplasia Using 3D MR Fingerprinting*.  
Ann Neurol, 2024, 96: 944-957.  
[https://doi.org/10.1002/ana.27049](https://doi.org/10.1002/ana.27049)

---

## Pipeline Overview
This pipeline outlines the preprocessing, feature generation, and analysis steps used in our study.

### 1. Data Organization
Organize your data in the following structure:

```
Data/
├── Patients/
│ ├── Patient01.nii
│ ├── Patient02.nii
│ └── ...
├── Target_subject/
│ └── Patient01.nii
├── HCs_test/
│ ├── V01.nii
│ ├── V02.nii
│ └── ...
├── DCs_test/
│ ├── V01.nii
│ ├── V02.nii
│ └── ...
├── HCs_norm/ # for normalization
│ ├── V01.nii
│ ├── V02.nii
│ └── ...

```

### 2. Image Registration
- Run `regis_2_MNI.sh` to register the MRF images into MNI space.

### 3. Tissue Segmentation
- Run `Fast_process_in_batch.sh` to segment the GM, WM, and CSF maps.

### 4. Normalization Data Generation
- Run `HCs_data_gen_4_norm.m` to generate the data for normalization.

### 5. Morphometric Maps (Optional)
- Run `MAP18` to create extension, thickness, and junction maps if needed.

### 6. Feature Extraction
- Run `MRF_vxl_2D_ft_gen.m` to normalize and create the 2D-level data for patients and controls (excluding lesional ROIs).  
- If running **FCD subtyping**:  
  - Set `subtype_flag = 1` to calculate additional features from MRF maps (entropy, uniformity).  
  - It is also recommended to set `VBM_flag = 1` to calculate VBM features during subgroup analysis.  

### 7. Statistical Analysis
- Run `Stats_ana.m` to perform statistical analysis between different groups.

### 8. Classification
- Run `Classification.m` to perform classification between two groups.  
- ⚠️ Please revise the group definitions according to your needs.  

---
