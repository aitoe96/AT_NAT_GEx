# AT_NAT_GEx
This repository contains the script eco_fun.R, used for analyzing whole-genome gene expression patterns in Arabidopsis thaliana accessions from an altitudinal gradient in southern Spain. The analysis focuses on understanding gene expression variations across diurnal, seasonal, and annual timescales under natural environmental conditions.

# Abstract

Understanding the genes, regulatory networks, and environmental cues involved in plant development is crucial in plant sciences. We conducted common garden experiments to quantify whole-genome gene expression patterns over diurnal, seasonal, and annual timescales across the life-cycle phenology of Arabidopsis thaliana with wild Iberian accessions. Key findings include:

    - Diurnal and seasonal comparisons structured gene expression patterns.
    - Accessions from contrasting environments exhibited variation in fitness-related traits.
    - Common differentially expressed genes and annotated genes across biological functions were identified.
    - Accessions with similar flowering times showed similar gene expression patterns, implicating plasticity in life-cycle phenology.
    - Harsh environmental conditions significantly impacted gene expression patterns.

# Repository Contents

   ## R Scripts:
        - eco_fun.R: The main script covering data loading, analysis, and result generation.

  ## Data Files:
        - ecofun_gene_expression.tsv: Gene expression data before normalization.
        - ecofun_gene_expression_df.tsv: Normalized gene expression data.
        - ecofun_design.tsv: Experimental design file.
   

# Setup Instructions

  ## Dependencies: 
  Ensure you have installed the necessary R packages listed in the script (ballgown, limma, clusterProfiler, etc.).

  ## File Paths: 
  Update input_path and output_path variables in the scripts to reflect your local directory structure.


# Results

## Generated results include:

    Differential gene expression analysis results in results/deg/ directory.
    Principal Component Analysis (PCA) results in results/pca/ directory.
    Various plots and figures related to the analysis.

# Conclusion

Our study highlights the plasticity in life-cycle phenology and gene expression exhibited by natural Arabidopsis thaliana accessions under the experimental conditions. 

# Contact

For questions or issues regarding this repository, please contact at aitoe96lg@gmail.com.



