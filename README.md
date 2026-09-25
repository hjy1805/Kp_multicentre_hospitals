# Repository of Genomic Signatures and Prediction of Clinical Severity in *Klebsiella pneumoniae* infections in a Multicenter Cohort

This repository contains all analysis code (bash, R, and Python) together with intermediate genomics and machine-learning artefacts generated for our multicentre large-scale genomic study of *Klebsiella pneumoniae*. Input datasets required to reproduce the machine-learning analyses are provided on the GitHub Releases page.

*Klebsiella pneumoniae* is a major causative agent of hospital-acquired infections worldwide, contributing substantially to morbidity, mortality, and healthcare burden. The clinical severity of K. pneumoniae infections is shaped by a complex interplay between antimicrobial resistance, virulence, and the underlying genomic background of infecting strains.  Despite extensive characterization of the genetic determinants of multidrug resistance and hypervirulence, the relationship between the genetic repertoire of K. pneumoniae and the clinical severity of infections remains inadequately understood.
## Methods
We analysed a large-scale nationwide collection of *K. pneumoniae* complex strains retrieved over seven years from five centres across the Kingdom of Saudi Arabia. Using comprehensive patient-level clinical data, we employed regression analyses, bacterial genome-wide association study (bGWAS) methods, and machine-learning approaches to elucidate the clinical significance of extended-spectrum β-lactamase/carbapenemase-producing (ESBL/CP) and hypervirulent *K. pneumoniae* (hvKp). hvKp was defined by the concurrent presence of all five hallmark biomarkers (*iucA*, *iroB*, *peg-344*, *rmpA*, and *rmpA2*), with isolates additionally carrying ESBL and/or CP determinants classified as convergent hvKp. We examined in-hospital all-cause mortality, intensive care unit (ICU) admission, and length of hospitalisation (LOS), identified genome-wide determinants associated with these outcomes, and employed machine-learning approaches to predict them using genomic biomarkers together with clinical metadata.
## Results
Infections caused by ESBL(+)/CP(+)-only strains showed consistently poorer outcomes than ESBL/CP-negative non-hvKp strains, including increased adjusted mortality hazard (HR 1.34, 95% CI 1.01–1.78), significantly poorer survival, and ~45% longer LOS. In contrast, hvKp-only infections showed no evidence of increased clinical severity, while only five isolates were classified as convergent hvKp. Carbapenem-resistance determinants were consistently associated with adverse outcomes after confounder adjustment, including higher odds of mortality (adjusted OR 2.20) and ICU admission (adjusted OR 1.75). bGWAS identified severity-associated accessory genes involved in carbohydrate metabolism, a type VI secretion system (T6SS) component, metabolic adaptation, and stress tolerance/persistence. Machine-learning models combining genomic and clinical predictors yielded average AUCs of 0.78 and 0.79 for mortality and ICU admission, respectively, and a moderate correlation between observed and predicted LOS (r = 0.59) in held-out test data.
## Conclusion
This study identifies key genomic determinants associated with severe K. pneumoniae infections, with antimicrobial-resistance markers showing the most consistent associations with adverse clinical outcomes. The predictive information provided by genomic biomarkers for mortality, ICU admission, and LOS highlights their potential value for improving clinical risk stratification and informing targeted infection-prevention strategies.

 



## Repository Structure

```plaintext
.
├── Bash_code                                     # Folder containing bash scripts
|   ├── add_date_city.sh                              # Bash script for adding the collection and location info to the alignment file
│   ├── annotate_prokka.sh                            # Assembled genome annotation bash script
│   ├── assembly_qc_quast.sh                          # Genome assembly quality assessment bash script
│   ├── assembly_unicycler.sh                         # Short-read genome assembly bash script
│   ├── assembly_unicycler_hybird.sh                  # Long-read hybrid genome assembly bash script
│   ├── call_snp_sites.sh                             # Bash script for calling SNPs from alignment by SNP-sites 
│   ├── mapping_snippy.sh                             # Bash script for mapping short reads against reference genome
│   ├── profiling_metaphlan.sh                        # Bash script for profiling the sequencing reads species as QC for contamination
│   ├── run_beast1.sh                                 # Bash script for running Bayesian Evolutionary Analysis Sampling Trees v1 (BEAST1)
│   ├── run_gubbins.sh                                # Bash script for filtering out polymorphic sites
│   ├── scan_amrfinder.sh                             # Bash script for AMR gene & virulence factor detection by AMRFinderPlus
│   └── typing_Kleborate.sh                           # Bash script for multi-function profiling of Klebsiella genome by Kleborate
├── R_code                                        # Folder containing R script
│   ├── global_samples.R                              # R script for clustering the external public isolates
│   ├── phylodynamic.R                                # R script for analysing phylodynamic results
│   ├── plasmid.R                                     # R script for clustering and visualising the plasmids 
│   └── seq_meta.R                                    #  R script for processing sequence metadata
├── Python_code                                     # Folder containing R script
│   ├── MLModels.ipynb                                # Python script for the machine learning models
│   ├── requirements.txt                              # Library requirements
│   ├── Mortality_XGBoost_results.csv                  # Results of the XGBoost model trained for Mortality label
│   ├── ICU_XGBoost_results.csv                        # Results of the XGBoost model trained for ICU Admission label 
│   └── LOS_NGBoost_results.csv                        # Results of the NGBoost model trained for LOS label 
├── Files                                         # Folder containing metadata and plasmid sequences
│   │ Plasmid_sequences                               # Folder containing plasmid sequences in fasta format
│   └──    └── [plasmid sequences in fasta format]
│   └──    └── GWAS Results for SNP and PanGenome
│   └──    └── peg-344 sequence and profile files 
└── README.md                             

```

##  Data Release structure
The large data release can be accessed through the release page or link: https://github.com/hjy1805/Kp_multicentre_hospitals/releases/tag/v1

```plaintext
DataSubmit/
├── ML                                              # Folder contains files for machine learning model training 
│   ├── Labels                                           # Folder of labels
│   │   ├──df_phenotype_ICU.csv                               # dataframe of ICU label of patients
│   │   ├──df_phenotype_LOS.csv                               # dataframe of Length of Stay (LOS) label of patients
│   │   └──df_phenotype_Mortality.csv                         # dataframe of mortality label of patients
│   ├── Predictors                                       # Folders of predictors
│   │   ├──df_phenotype_ICU.csv                               # dataframe of ICU predictors of patients
│   │   ├──df_phenotype_LOS.csv                               # dataframe of Length of Stay (LOS) predictors of patients
│   │   └──df_phenotype_Mortality.csv                         # dataframe of mortality predictors of patients
├── Plasmid                                         # Folder contains files for plasmids
│   ├── plasmid_profiles.csv                             # plasmid profiles for AMR, virulence genes and replicon
│   └── PlasmidONTAccession.csv                          # ENA and GenBank accessions for ONT-sequenced plasmid origin samples
├── gene_presence_absence_PanGenome.csv             # Panaroo gene presence–absence matrix used for pangenome analysis in csv format 
└── pan_genome_reference.fa                         # Panaroo pangenome reference sequences in FASTA format    

```


## Reference

If you use this repository, please cite the corresponding publication (details to be added upon acceptance).

---

## Contact

For questions, issues, or collaboration related to this repository, please contact:

**Jiayi Huang, PhD Student**  
Infectious Disease Epidemiology Laboratory  
Biomedical Sciences Division (BioMed) 
King Abdullah University of Science and Technology (KAUST)  

📧 **Email:** jiayi.huang@kaust.edu.sa

**Sara Iftikhar, Research Assistant**  
Infectious Disease Epidemiology Laboratory  
Biomedical Sciences Division (BioMed)
King Abdullah University of Science and Technology (KAUST)  

📧 **Email:** sara.iftikhar@kaust.edu.sa
