# Genomic Signatures and Prediction of Clinical Severity in *Klebsiella pneumoniae* infections in a Multicenter Cohort


## Abstract

*Klebsiella pneumoniae* is a major causative agent of hospital-acquired infections. The emergence of strains that combine resistance to last-resort antimicrobials with hypervirulence has become a pressing public-health challenge. Despite extensive characterisation of the genetic determinants of multidrug resistance and hypervirulence, the relationship between the genetic repertoire of *K. pneumoniae* and the clinical severity of infections remains inadequately understood.

## Methods
We analysed a nationwide large-scale collection of 1,319 *K. pneumoniae* complex strains retrieved over eight years from five centres across the Kingdom of Saudi Arabia. Using detailed and comprehensive clinical metadata, we employed a range of regression analyses, genome-wide association study methods and machine-learning techniques to decipher the clinical significance of ESBL/carbapenemase-producing positive (ESBL/CP), hypervirulent (hv), and convergent resistant and hv (ESBL/CP+/hv) *K. pneumoniae* strains. We examined clinical severity outcomes, including in-hospital mortality rate, ICU admission rate and length of hospitalisation (LOS) across these *K. pneumoniae* types, identified genome-wide determinants linked with clinical severity and used genomic biomarkers together with clinical metadata to predict clinical severity outcomes.

## Results
Infections caused by convergent strains exhibited the highest severity, with nearly double the in-hospital mortality (reaching 42% at 90 days), a 2.1-fold greater likelihood of ICU admission (p < 0.001), and an average of 150 more hospitalisation days compared to the other infections—indicating an additive effect of hypervirulence and multidrug resistance. Carbapenem resistance determinants showed the strongest association with adverse outcomes, even after adjusting for the presence of other resistance and virulence genes and clinical confounder features. The GWAS analysis revealed associations of the clinical outcomes with accessory genes involved in carbohydrate metabolism and the Type VI secretion system (T6SS) machinery, metabolic-adaptation and stress-tolerance/persistence loci. Additional significant associations were identified with SNPs in ABC-transporters, cell-envelope systems, sugar transporter families and RND-family efflux systems. Machine-learning models yielded Area Under the Curve (AUC) values of 0.80 and 0.87 for mortality and ICU admission, respectively, and explained 65% of the variance in LOS on unseen data when trained using genomic information—showing that including genomic data significantly improved predictions over clinical data alone.

## Conclusion
These findings provide quantitative measures of genetic risk factors for *K. pneumoniae* infections and highlight the role of carbapenem resistance as the main driver of clinical severity. They also underscore the potential of genomic biomarkers as predictive diagnostic tools for clinical management and infection-prevention strategies for *K. pneumoniae*.

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
├── Python_code                                        # Folder containing R script
│   ├── MLModels.ipynb                                # Python script for the machine learning models
│   ├── requirements.txt                                # Library requirements
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


## Contacts
For inquiries regarding this research, please contact:

Jiayi Huang

Email: jiayi.huang@kaust.edu.sa

PhD student

Infectious Disease Epidemiology Lab

Biological and Environmental Science and Engineering (BESE) Division

King Abdullah University of Science and Technology (KAUST)
