# Genomic Signatures and Prediction of Clinical Severity in *Klebsiella pneumoniae* Infections in a Multicentre Cohort

Analysis code and supporting genomic and machine-learning artefacts for our multicentre study of *K. pneumoniae*. The manuscript is under revision following peer review.

[Preprint](https://doi.org/10.64898/2026.02.02.26345332) · [Data release](https://github.com/hjy1805/Kp_multicentre_hospitals/releases/tag/v1) · [GenoPredict](https://genopredict.kaust.edu.sa)

[Study overview](#study-overview) · [Requirements](#requirements) · [Running the analyses](#running-the-analyses) · [Repository structure](#repository-structure) · [Citation](#citation)

## Study Overview

We investigate how antimicrobial resistance, virulence, and genomic background relate to the clinical severity of *K. pneumoniae* infections.

- **Cohort:** a nationwide collection of *K. pneumoniae* complex strains collected over seven years at five centres in Saudi Arabia, with patient-level clinical data.
- **Outcomes:** in-hospital all-cause mortality, intensive care unit (ICU) admission, and length of hospitalisation (LOS).
- **Analyses:** regression, bacterial genome-wide association studies (bGWAS), and machine learning using genomic biomarkers and clinical metadata.
- **Pathotypes:** extended-spectrum β-lactamase/carbapenemase-producing (ESBL/CP), hypervirulent (hvKp), and convergent hvKp strains. hvKp requires all five biomarkers (*iucA*, *iroB*, *peg-344*, *rmpA*, and *rmpA2*); convergent hvKp additionally carries ESBL and/or CP determinants.

### Key Findings

- **ESBL(+)/CP(+)-only infections:** increased adjusted mortality hazard (HR 1.34, 95% CI 1.01–1.78), significantly poorer survival, and approximately 45% longer LOS than ESBL/CP-negative non-hvKp infections.
- **hvKp-only infections:** no evidence of increased clinical severity; only five isolates were classified as convergent hvKp.
- **Carbapenem-resistance determinants:** associated with higher odds of mortality (adjusted OR 2.20) and ICU admission (adjusted OR 1.75) after confounder adjustment.
- **bGWAS:** severity-associated accessory genes involved in carbohydrate metabolism, a type VI secretion system component, metabolic adaptation, and stress tolerance/persistence.

Models combining genomic and clinical predictors achieved the following performance in held-out test data:

| Outcome | Performance |
|---|---|
| Mortality | Average AUC: 0.78 |
| ICU admission | Average AUC: 0.79 |
| Length of stay | Observed–predicted correlation: *r* = 0.59 |

Antimicrobial-resistance markers showed the most consistent associations with adverse outcomes, highlighting the potential of genomic biomarkers for clinical risk stratification and targeted infection prevention.

## GenoPredict

[GenoPredict](https://genopredict.kaust.edu.sa) is our portal for clinical severity prediction from *K. pneumoniae* genomes.

1. **Upload** an assembled bacterial genome.
2. **Scan** for genomic biomarkers: unitigs and antimicrobial resistance (AMR) and virulence genes.
3. **Receive** estimates of mortality risk, ICU admission risk, and length of stay, each with an approximate 95% interval.

[![GenoPredict portal for clinical severity prediction from K. pneumoniae genomes](docs/images/genopredict.png.jpeg)](https://genopredict.kaust.edu.sa)

## Repository Structure

| Directory | Contents |
|---|---|
| [bash_code/](bash_code/) | Genome assembly, quality control, annotation, genomic profiling, and phylogenetic analysis scripts |
| [R_code/](R_code/) | Statistical analyses, transmission analyses, metadata processing, and visualisation |
| [python_code/](python_code/) | Machine-learning notebook, dependencies, and model results |
| [Files/](Files/) | Sample metadata, plasmid sequences, GWAS results, and *peg-344* files |
| [docs/images/](docs/images/) | GenoPredict portal screenshot used in this README |

## Requirements

### Python

- **Python:** the notebook records Python 3.11.7 in its metadata.
- **Scientific computing:** NumPy, pandas, and SciPy.
- **Modelling:** scikit-learn, XGBoost, and NGBoost.
- **Interpretation and plotting:** SHAP, Matplotlib, and seaborn.
- **Utilities:** joblib, tqdm, numba, and JupyterLab.

Package constraints are provided in [python_code/requirements.txt](python_code/requirements.txt). From the repository root, create an environment and install them:

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -r python_code/requirements.txt
```

### R

Use R with the packages required by the selected analysis.

| Analysis | Packages loaded by the scripts |
|---|---|
| Clinical outcomes and model evaluation | `readr`, `dplyr`, `tidyr`, `stringr`, `forcats`, `ggplot2`, `purrr`, `scales`, `survival`, `broom`, `forestmodel`, `cmprsk`, `survminer` |
| Transmission and threshold analyses | `ape`, `dplyr`, `tidyr`, `purrr`, `tibble`, `readr`, `stringr`, `ggplot2`, `adegenet`, `igraph`, `gt` |
| Metadata and plasmid analyses | `tidyverse`, `readxl`, `ggrepel`, `factoextra`, `fpc`, `RColorBrewer`, `Rtsne` |
| Global isolate analyses | `tidyverse`, `ape`, `treedataverse`, `ggnewscale`, `rhierbaps`, `phytools`, `ggsci`, `aplot`, `ggimage`, `rsvg`, `cutpointr`, `ggsignif` |

For the clinical outcome and model-evaluation scripts, install the following from an R console:

```r
install.packages(c(
  "readr", "dplyr", "tidyr", "stringr", "forcats", "ggplot2",
  "purrr", "scales", "survival", "broom", "forestmodel", "cmprsk", "survminer"
))
```

## Running the Analyses

### 1. Prepare the repository and data

```bash
git clone https://github.com/hjy1805/Kp_multicentre_hospitals.git
cd Kp_multicentre_hospitals
```

- Download and extract the [data release](https://github.com/hjy1805/Kp_multicentre_hospitals/releases/tag/v1).
- Install the Python or R dependencies above for the analysis you plan to run.
- Run R scripts from the repository root.

### 2. Python: machine-learning analyses

With the Python environment activated, launch the notebook from the repository root:

```bash
python -m jupyterlab python_code/MLModels.ipynb
```

The notebook contains sections for logistic regression, XGBoost, quantile regression, NGBoost, and SHAP interpretation.

### 3. R: statistical and genomic analyses

| Analysis | Script |
|---|---|
| Mortality with discharge as a competing event | [competing_risk_pathotype.R](R_code/competing_risk_pathotype.R) |
| Cox regression by pathotype | [cox_pathotype_mortality.R](R_code/cox_pathotype_mortality.R) |
| Mortality and ICU determinants | [mortality_determinants.R](R_code/mortality_determinants.R), [icu_determinants.R](R_code/icu_determinants.R) |
| Transmission networks | [transmission_network.R](R_code/transmission_network.R) |
| K-mer threshold and sensitivity analyses | [kmer_transmission_analysis.R](R_code/kmer_transmission_analysis.R), [kmer_threshold_analysis.R](R_code/kmer_threshold_analysis.R) |

### 4. Bash: genome-processing workflows

The Bash scripts include cluster-specific paths, software modules, and Slurm settings. Install or load the tools used by the selected script and adjust its paths and resource requests before submission.

For example, `assembly_qc_quast.sh` requires QUAST, a `tag_list` file, and assemblies at `<sample>/hybrid_assembly/assembly.fasta`. From the directory containing those inputs, submit the configured script on a Slurm cluster:

```bash
sbatch /path/to/Kp_multicentre_hospitals/bash_code/assembly_qc_quast.sh
```

Set the job-array range to the number of samples in `tag_list`.

## Data Release Structure

Datasets are available in the [GitHub data release (v1)](https://github.com/hjy1805/Kp_multicentre_hospitals/releases/tag/v1). Paths below are relative to `DataSubmit/`.

| Path | Contents |
|---|---|
| `ML/Labels/` | Outcome labels for mortality, ICU admission, and length of stay |
| `ML/Predictors/` | Predictor datasets for each outcome |
| `Plasmid/plasmid_profiles.csv` | Plasmid profiles for AMR genes, virulence genes, and replicons |
| `Plasmid/PlasmidONTAccession.csv` | ENA and GenBank accessions for the source samples of ONT-sequenced plasmids |
| `gene_presence_absence_PanGenome.csv` | Panaroo pangenome gene presence–absence matrix |
| `pan_genome_reference.fa` | Panaroo pangenome reference sequences in FASTA format |

## Citation

If you use this repository in your work, please cite our preprint:

```bibtex
@article{Malaikah2026.02.02.26345332,
  author = {Malaikah, Mohammad and Alyami, Rfeef Yousf and Huang, Jiayi and Fallatah, Omniya and Milner, Mathew and Zhou, Ge and Hirayban, Raneem and Iftikhar, Sara and Banzhaf, Manuel and Li, Yan and Senok, Abiola and Hala, Sharif Matouq and Batook, Nouf and Alsharif, Dhuha and AlShahrani, Alhanouf Sultan and Alamri, Abdulfattah Wasel and AlJohani, Sameera M and Kaaki, Mai M and Alwan, Basaam and Absar, Mohammad and Ali, Mohammad Elsir Al. and Sadah, Haitham S and Zakri, Samer and Bosaeed, Mohammad and Pain, Arnab and Moradigaravand, Danesh},
  title = {Genomic Signatures and Prediction of Clinical Severity in Klebsiella pneumoniae infections in a Multicenter Cohort},
  journal = {medRxiv},
  year = {2026},
  doi = {10.64898/2026.02.02.26345332},
  url = {https://www.medrxiv.org/content/early/2026/02/03/2026.02.02.26345332}
}
```


## Contact

For questions, issues, or collaboration:

| Contact | Role | Email |
|---|---|---|
| Jiayi Huang | PhD Student | [jiayi.huang@kaust.edu.sa](mailto:jiayi.huang@kaust.edu.sa) |
| Sara Iftikhar | Research Assistant | [sara.iftikhar@kaust.edu.sa](mailto:sara.iftikhar@kaust.edu.sa) |

Infectious Disease Epidemiology Laboratory, Biomedical Sciences Division (BioMed), King Abdullah University of Science and Technology (KAUST).
