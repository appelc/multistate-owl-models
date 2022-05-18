# multistate-owl-models
Multi-state occupancy models for northern spotted owls

-------------------------

## FOLDER STRUCTURE

### bugs
text files with models for JAGS (see model descriptions, below)

### data2018
- nsod_dh_2018.csv: spotted owl detection history
- covariates folder: contains raw and standardized versions of covariate data

### scripts
- step01_models: folder with R scripts to run each model (see descriptions, below)
- step02_submodel_selection.R: evaluating sub-models using sequential-by-submodel selection
- step03_model_evaluation.R: evaluating most-supported model from sub-model selection (includes derived parameters for reporting in-text and tables, creating marginal plots, etc.)
- step04_posthoc_model.R: evaluating post-hoc model for interpreting covariate effects on psi and R