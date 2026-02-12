# SSRI Organoids – Single-Cell Analysis of Sertraline Exposure

**Research Article:** (add link to paper once published)  
**Project:** Single-cell RNA sequencing analysis of human brain organoids exposed to sertraline (SSRI)

---

## 🧪 Overview
This repository contains analysis workflows and scripts used to study the effects of **sertraline (a selective serotonin reuptake inhibitor, SSRI)** on human brain organoids using single-cell RNA sequencing (scRNA-seq).

---

## 📁 Repository Structure
1_QC_dimred_clustering.R
- Quality control 
- Dimensionality reduction 
- Cell type annotation

2_DE_analysis.R
- Pseudobulk Differential gene expression

3_OrganoidAtlas_ref_mapping.ipynb
- Reference mapping to an external neural organoid  atlas (HNOCA)

5_composition.R
- Cell type composition analysis using cacoa.

8_trajectory.R
- `slingshot` for lineage inference
- Pseudotime estimation
- Visualization of lineage structures in UMAP space

6_disorder_enrich.R
- Overlap analysis with curated disorder gene sets

Contact
Sellgren Lab (carl.sellgren@ki.se)
Applied Developmental Neurobiology
Karolinska Institute, Sweden

For questions regarding this repository, please open an issue.

