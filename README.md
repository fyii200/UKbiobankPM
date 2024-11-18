### Source code is in the directory named 'code'

Order in which the scripts were executed: 
1. **randomSampler.py** : Python script for randomly sampling 2,000 participants in the UK Biobank with high myopia (spherical equivalent refraction ≤ -5D) in at least one eye.
2. **preprocessData.R** : R script for preprocessing relevant variables (explored risk factors, fundus imaging features & pathologic myopia labels) analysed in the study.
3. **analysis.R**       : R script for performing the main analysis (prevalence estimates and logistic regression).

Note 1: **utils.R** is a non-executable R script containing various functions used in 'analysis.R'.

Note 2: Source code for deriving fundus imaging features is available [elsewhere](https://github.com/fyii200/MyopiaRetinalFeatures).

