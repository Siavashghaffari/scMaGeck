# scMaGeck
This repo contains:
1. A python notebook to prepare all the tools for eff_scmageck run: `scMaGeck_EFF_prep.ipynb`
2. Data conversion file from 10x to Seurat in R to be a proper file to input eff_scmageck run: `data_conversion.R`
3. eff_scmageck run on all genes gene by gene: `EFF_run.R`
4. eff_scmageck run on a list of genes: `EFF_GeneList.R`
5. eff_scmageck run on some target gens of interest: `EFF_targetGenes.R`

- First execute steps 1 and 2. Following that proceed with either 3, 4 or 5 based on specific task and objective you wish to accomplish

## input folder
Should contain:
1. the main annData which contains the raw counts
2. the guide RNA AnnData
3. A csv file which has information about gRNA sequences (the file that we used for cumulus is a good example)

## Outputs
1. the PS score matrix containing the PS scores of each cells for each perturbed gene 
2. A list of lr results, including: beta score (as a data.frame), the associated p value (as a data.frame), and the regression matrix that is used for linear regression to test the association of each perturbed gene with all possible genes


## Authors and acknowledgment
This work was developed by Siavash Ghaffari. For any questions, feedback, or additional information, please feel free to reach out. Your input is highly valued and will help improve and refine this pipeline further.