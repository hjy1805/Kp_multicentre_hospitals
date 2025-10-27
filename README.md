# Unravelling the epidemiological landscape of Klebsiella pneumoniae in multicentre hospitals

## Abstract

Antimicrobial resistance (AMR) in Klebsiella pneumoniae has become a major healthcare threat, particularly within tertiary-care hospitals in Saudi Arabia where carbapenem-resistant and hypervirulent lineages increasingly converge. This multi-center, multi-omics investigation, conducted across five tertiary-care hospitals between 2018 and 2024, integrated antimicrobial susceptibility testing, short- and long-read whole-genome sequencing, plasmidome reconstruction, and lineage-focused RNA sequencing with patient metadata and clinical outcomes, including survival and multivariable Cox regression analyses. Among 1,319 sequenced isolates, we observed high allelic diversity (404 sequence types) but oligoclonal dominance of high-risk lineages—most notably ST2096, ST14, and ST147—which displayed distinct regional distributions. Convergent multidrug-resistant–hypervirulent (MDR-hv) infections rose from <5% to nearly 20% over the study period and were independently associated with prolonged hospitalization, greater ICU admission, and higher 30-, 60-, and 90-day mortality. Long-read assemblies (n = 179) revealed a plasmidome dominated by IncFIB(K), IncFII(K), and Col-type replicons, with frequent ColKP3–IncF co-carriage and occasional large IncHI1B–IncF hybrids encoding both resistance and virulence determinants. Canonical backbone–enzyme associations—IncL/M–blaOXA-48-like, IncC/IncFII(K)–blaNDM, and ColKP3–blaOXA-232—illustrated modular gene flow among epidemic lineages. Comparative transcriptomics of 43 representative isolates identified distinct lineage-specific expression programs: ST14exhibited activation of anaerobic respiration and organic-acid catabolism pathways; ST147 up-regulated oxidative and small-molecule metabolic networks; and ST2096 enriched histidine, imidazole, and xenobiotic-response modules, suggesting complementary adaptive strategies that may facilitate persistence, stress tolerance, and plasmid maintenance in clinical niches. Collectively, these findings demonstrate that K. pneumoniae risk in Saudi tertiary-care hospitals is increasingly driven by convergent MDR-hv lineages—particularly ST2096—with measurable clinical consequences. The integration of plasmid structure, lineage-specific expression profiles, and patient outcomes supports a biologically coherent model of convergent evolution and underscores the urgent need for targeted antimicrobial stewardship, genomic surveillance audits, and infection-prevention strategies tailored to Saudi healthcare systems.

## File Structure

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
├── R&Python_code                                 # Folder containing R script and Jupyter notebook
│   ├── metadata_processing.R                         # R script for the relevant metadata processing
│   └── plasmid_pygenomeviz.ipynb                     # Python script for the plasmid alignment and visualization by pyGenomeViz
├── Files                                         # Folder containing metadata and annotated plasmid sequences
│   ├── Annotated_plasmid                             # Folder containing annotated plasmid sequences in gbk format
│   │   └── [Annotated plasmid files in gbk format]   
│   └── Samples_metadata.csv                          # Metadata of samples that were in-house sequenced in this study
└── README.md                             

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
