# Repository of Genomic Signatures and Prediction of Clinical Severity in *Klebsiella pneumoniae* infections in a Multicenter Cohort

This repository contains all analysis code (bash, R, and Python) together with intermediate genomics and machine-learning artefacts generated for our multicentre large-scale genomic study of *Klebsiella pneumoniae*. Input datasets required to reproduce the machine-learning analyses are provided on the GitHub Releases page.

## Abstract

*Klebsiella pneumoniae* is a major causative agent of hospital-acquired infections worldwide, contributing substantially to morbidity, mortality, and healthcare burden.. The emergence of strains that combine resistance to last-resort antimicrobials with hypervirulence has become a pressing public-health challenge. Despite extensive characterization of the genetic determinants of multidrug resistance and hypervirulence, the relationship between the genetic repertoire of *K. pneumoniae* and the clinical severity of infections remains inadequately understood.

## Methods
We analyzed a nationwide large-scale collection of 1,306 *K. pneumoniae* complex strains retrieved over seven years from five centres across the Kingdom of Saudi Arabia. Using detailed and comprehensive patient-level clinical data, We employed a range of regression analyses, genome-wide association study (GWAS) methods, and machine-learning approaches to elucidate the clinical significance of ESBL/carbapenemase-producing (ESBL/CP), hypervirulent, and convergent ESBL(+)/CP(+) hypervirulent K. pneumoniae strains. We examined clinical severity outcomes including in-hospital all-cause mortality rate, ICU admission rate and length of hospitalisation (LOS) across these K. pneumoniae types, identified genome-wide determinants linked with clinical severity and used machine learning approaches to predict clinical severity outcomes from genomic biomarkers together with clinical metadata.

## Results
Infections caused by convergent strains exhibited the greatest clinical severity, showing nearly double the in-hospital mortality (reaching 42% at 90 days), a 2.4-fold higher likelihood of ICU admission, and an average 150% increase in LOS compared to infections caused by susceptible and non-hypervirulent strains. Our findings indicate an additive effect of hypervirulence and multidrug resistance on disease severity. Carbapenem resistance determinants showed the strongest association with adverse outcomes, even after adjusting for the presentce of other resistance and virulence genes and clinical confounder features. The GWAS analysis revealed associations of the clinical outcomes with accessory genes involved in carbohydrate metabolism and the Type VI secretion system (T6SS) machinery, metabolic-adaptation and stress-tolerance/persistence loci. Additional significant associations were identified with SNPs in ABC-transporters, cell-envelope systems, sugar transporter families and RND-family efflux systems. Machine-learning models yielded average Area Under the Curve (AUC) values of 0.78 and 0.79 for mortality and ICU admission, respectively, and exhibited strong monotonic association between observed and predicted outcomes for LOS, with an average correlation of 0.59 on unseen test data when trained using combined genomic and clinical predictors.

## Conclusion
This study identifies key genomic determinants that drive severe *K. pneumoniae* infections, with carbapenem-resistance markers emerging as the leading contributors to poor clinical outcomes. The strong predictive performance of genomic biomarkers, particularly for mortality, ICU admission, and LOS, highlights their value in enhancing diagnostic precision, improving clinical risk stratification, and informing targeted infection-prevention strategies. 



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
└── README.md                             

```

##  Data Release structure
The large data release can be accessed through the release page or link: https://github.com/hjy1805/Kp_multicentre_hospitals/releases/tag/v1

```plaintext
DataSubmit/
├── GWAS                                             # Folder GWAS results
|   ├── Death_GWAS_PanGenome.csv                         # GWAS results for patient mortality using a pangenome presence–absence matrix
│   ├── ICU_GWAS_PanGenome.csv                           # GWAS results for patient ICU using a pangenome presence–absence matrix
│   ├── LOS_GWAS_PanGenome.csv                           # GWAS results for patient Length of Stay (LOS) using a pangenome presence–absence matrix
│   ├── SNPs_ICU_GWAS.csv                                # GWAS results for patient ICU using a single-nucleotide polymorphisms (SNPs)
│   ├── SNPs_LOS_GWAS.csv                                # GWAS results for patient Length of Stay (LOS) using a single-nucleotide polymorphisms (SNPs)
│   └── SNPs_Mortality_GWAS.csv                          # GWAS results for patient mortality using a single-nucleotide polymorphisms (SNPs)
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
Biological and Environmental Science and Engineering (BESE) Division  
King Abdullah University of Science and Technology (KAUST)  

📧 **Email:** jiayi.huang@kaust.edu.sa
