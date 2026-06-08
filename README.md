# Genome-wide Scans for  Locus Association with Antimalarial Drug Responses in *Plasmodium falciparum*

This repository provides scripts to perform **genome-wide scans for signatures of selection and recombination** in *Plasmodium falciparum* using population genomic data.
It also includes scripts for **mixed effets-based genotype-phenotype association analysis**, the **temporal trend analysis of haplotype prevalence** and **haplotype frequency plotting**. 
The workflows implemented here were developed in the context of the **PX1 manuscript**, but they are generalizable and can be applied to other datasets with minimal adaptation.


---

## Features

This repository contains **two main genomic analysis pipelines**, implemented in two files: [`Selection.R`](https://github.com/Karaniare/Genome_wide_scan_for_Selection_and_Recombination_PX1/blob/main/PX1_Manuscript_github.R) and [`making_ldhat.sh`](https://github.com/Karaniare/Genome_wide_scan_for_Selection_and_Recombination_PX1/blob/main/making_ldhat.sh).  
A genotype-phenotype association analysis pipeline, structured as a single file [`MIxed-Effects_Model.R`](https://github.com/Karaniare/PfWGS-DRscanner/blob/main/main/MIxed-Effects_Model.R) and a haplotype trend estimation pipeline file [`Temporal_Haplotype_Frequency_Analysis.R`](https://github.com/Karaniare/Genome_wide_scan_for_Selection_and_Recombination_PX1/blob/main/Temporal_Haplotype_Frequency_Analysis.R) 

Each file includes **step-by-step instructions** to guide users through the analyses. Open them to access the instructions and the codes.

---

### 1. Genome-wide Selection Scans  
**`PX1_Manuscript_github.R`**

This R script performs **two independent and complementary tests for positive selection**:

#### IBD-based selection scan
- Computes *iR statistics* using  
  👉 [IsoRelate](https://github.com/bahlolab/isoRelate)
- Based on excess identity-by-descent (IBD) sharing across the genome

#### Haplotype-based selection scan
- Computes **Integrated Haplotype Score (iHS)**
- Based on extended haplotype homozygosity (EHH)
- Implemented using the  
  👉 [`rehh`](https://cran.r-project.org/web/packages/rehh/index.html) R package

The script covers the full workflow, from input processing to genome-wide visualization of selection signals.

---

### 2. Genome-wide Recombination Rate Estimation  
**`making_ldhat.sh`**

This bash script estimates **recombination rates across chromosomes** using **LDhat**:

- Chromosome-wise processing
- Estimation of population-scaled recombination rates (ρ)
- Designed for large-scale, genome-wide analyses

Inline comments provide clear guidance on how to run and adapt the pipeline.

---

### 3. Genotype-Phenotype Assocation (Mixed-Effects Model)  
**`MIxed-Effects_Model.R`**

This R script test genotype-phenotype association based on mixed-effects model. It calculates:

- Cohen's effect size
- Conditional and Marginal R<sup>2<sup>
- Incremental R<sup>2<sup>
- Interaction terms between genes (P-values and incremental R<sup>2<sup>)

---

### 4. Haplotype Trend Analysis (Mann-Kendall Test)  
**`Temporal_Haplotype_Frequency_Analysis.R`**

This R script performs Mann-Kendall Test for monotonic trend and for plotting haplotype prevalences.

---

## Input Data

- For genomic analyses: a **genome-wide VCF file**
  - The VCF should be filtered for:
    - Genotype missingness
    - Minor allele frequency (MAF)
 
- For genotype-phenotype assocation: a simple **spreadsheet** containing at least:
  - PX1 haplotype,
  - K13 haplotype
  - Year of sample collection
  - Site
  - IC<sub>50<sub>
  - Drug assay name
  - Complexity of infection
  - Parasitemia
  - Within-sample allele frequency (optional, good for filtering)
 
- For Haplotype Trend Analysis: a simple **spreadsheet** containing at least:
  - PX1 haplotype,
  - K13 haplotype
  - Year of sample collection
  - Site
  - Within-sample allele frequency (optional, good for filtering) 

---

## Installation

Clone the repository:

```bash
git clone https://github.com/Karaniare/Genome_wide_scan_for_Selection_and_Recombination_PX1.git
cd Genome_wide_scan_for_Selection_and_Recombination_PX1```
```
## Resources

- **Example data for selection analysis**  
  *(https://github.com/Karaniare/Genome_wide_scan_for_Selection_and_Recombination_PX1/tree/main/Examples%20(selection))*

- **Example data for recombination rate estimation**  
  *(https://github.com/Karaniare/Genome_wide_scan_for_Selection_and_Recombination_PX1/tree/main/Examples)*

- **LDhat documentation and resources**
  - GitHub repository: https://github.com/auton1/LDhat
  - Includes source code, likelihood tables, and usage instructions

---

## Dependencies

### R Packages

The following R packages are required for the selection analyses:

- `isoRelate`
- `moimix`
- `vcfR`
- `ggplot2`
- `ggsci`
- `ggpubr`
- `rehh`
- `dplyr`
- `SeqArray`
- `R.utils`

### External Software

- **LDhat**  
  Required for genome-wide recombination rate estimation  
  - https://github.com/auton1/LDhat

---

## Notes

- The pipelines were developed for *Plasmodium falciparum* but can be adapted to other organisms.
- Users are encouraged to inspect and adjust filtering thresholds depending on sample size and study design.

---

## Citation

If you use this code, please cite the *([Niare et al. bioRxiv 2025, PMID: 40766679](https://pubmed.ncbi.nlm.nih.gov/40766679/))* and the relevant software packages (**IsoRelate**, **rehh**, and **LDhat**).



