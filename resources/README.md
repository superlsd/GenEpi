# Resources

This directory contains resources generated during this study, including harmonized disease nomenclature, gene and phenotype associations, molecular and epidemiological disease–disease networks, sex-specific disease trajectories, and drug–disease predictions. These files may support reproducibility, secondary analyses, and future studies.

## Contents

### A. `disease_name_consensus_official_names_PubMed.csv`

This table provides the matched disease name for each ICD-10 code, derived through fuzzy string matching and word-embedding approaches applied to multiple data sources.

- **Disease Name**: all identified alternative disease names and subdiagnoses.
- **Consensus Name**: the harmonized disease name used throughout the manuscript, obtained using a revised centroid-based approach.
- **PubMed_query_name**: a further-cleaned version of the consensus name used to query PubMed. The table also reports the number of retrieved publications for each disease.

### B. `allsources_icd_gene_dict_unified_v2026_official.pickle`

This pickle file contains a Python dictionary in which keys are ICD-10 codes and values are lists of associated genes (gene–disease associations, GDAs). Associations were harmonized across eight major databases using the disease nomenclature described in file A.

### C. `ICD_HPO_Genetic_dict_2026_official.pickle`

This pickle file contains a Python dictionary in which keys are ICD-10 codes and values are lists of associated Human Phenotype Ontology (HPO) terms. The terms were derived from gene enrichment between the associations in file B and genes reported in the Human Phenotype Ontology.

### D. `G_gene_icd_disp_significant_post_df_v2026_official.tsv`

This table contains edges (disease pairs) in the molecular disease–disease network (DDN), together with the strength of their molecular association (`weight`). Disease pairs represent statistically significant molecular overlap according to a Fisher’s exact test.

### E. `G_elma_icd_male_df_v2026_official.tsv`

This table contains disease pairs in the male-specific epidemiological disease–disease network (eDDN). Each pair was reported together in the Austrian population more frequently than expected.

### F. `G_elma_icd_female_df_v2026_official.tsv`

This table contains disease pairs in the female-specific epidemiological disease–disease network (eDDN). Each pair was reported together in the Austrian population more frequently than expected.

### G. `male_trajectories_dis_sort_bytime_grouped_all_add_evolving_dict_v2026_official.pickle`

This pickle file contains a Python dictionary of male-specific disease trajectories. Dictionary keys are trajectory identifiers, and values are disease lists ordered by chronological progression.

### H. `female_trajectories_dis_sort_bytime_grouped_all_add_evolving_dict_v2026_official.pickle`

This pickle file contains a Python dictionary of female-specific disease trajectories. Dictionary keys are trajectory identifiers, and values are disease lists ordered by chronological progression.

### I. `Supp_Tab5_disease_molecular_drugs_df_v2026_official_Pubmed.csv`

This table reports drug–disease predictions based on molecular overlap, together with the related PubMed publication counts reported in the literature.

### J. `SupplementaryTable11_ICD10Codes_Drug_TxAgent.xlsx`

This Excel workbook documents the prompting strategy and results for support of sex-specific disease trajectories and drug predictions.
