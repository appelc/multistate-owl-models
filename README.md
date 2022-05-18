# multistate-owl-models
Multi-state occupancy models for northern spotted owls

-------------------------

## CONTENTS

### bugs
- contains text files with models for JAGS (see model descriptions, below)

### data2018
- nsod_dh_2018.csv: spotted owl detection history
- covariates folder: contains raw and standardized versions of covariate data

### scripts
- step01_models: folder with R scripts to run each model (see descriptions, below)
- step02_submodel_selection.R: evaluating sub-models using sequential-by-submodel selection
- step03_model_evaluation.R: evaluating most-supported model from sub-model selection (includes derived parameters for reporting in-text and tables, creating marginal plots, etc.)
- step04_posthoc_model.R: evaluating post-hoc model for interpreting covariate effects on psi and R


## MODEL DESCRIPTIONS

01a: p submodel (initial)		  p(NOISE + EFFORT + SEASON + BO_WEEKLY + NR200) delta(.) psi(.) R(.)  
01b: p submodel (slab-spike)	p(NOISE + EFFORT + SEASON + BO_WEEKLY + NR200) delta(.) psi(.) R(.)

02a: psi submodel (initial)	    p(NOISE + EFFORT) delta(.) psi(STUDY_AREA + NR500 + BO_TOTAL) R(.)  
02b: psi submodel (slab-spike)	p(NOISE + EFFORT) delta(.) psi(STUDY_AREA + NR500 + BO_TOTAL) R(.)

03a: R submodel (initial)		  p(NOISE + EFFORT) delta(.) psi(.) R(STUDY_AREA + NR500 + BO_TOTAL)  
03b: R submodel (slab-spike)	p(NOISE + EFFORT) delta(.) psi(.) R(STUDY_AREA + NR500 + BO_TOTAL)

04: most-supported model from sub-model selection process
	p(NOISE + EFFORT) delta(.) psi(.) R(.)

05a: post-hoc model, "normal1" priors: dnorm(0,0.368)  
	p(NOISE + EFFORT) delta(.) psi(STUDY_AREA + NR500 + BO_TOTAL) R(STUDY_AREA + NR500 + BO_TOTAL)
	
05b: post-hoc model, "normal2" priors: dnorm(0,0.5)  
	p(NOISE + EFFORT) delta(.) psi(STUDY_AREA + NR500 + BO_TOTAL) R(STUDY_AREA + NR500 + BO_TOTAL)

05c: post-hoc model, "uniform" priors: dunif(-5,5)  
	p(NOISE + EFFORT) delta(.) psi(STUDY_AREA + NR500 + BO_TOTAL) R(STUDY_AREA + NR500 + BO_TOTAL)
	
06: post-hoc model, "normal1" priors: dnorm(0,0.368)  
	p(NOISE + EFFORT) delta(.) psi(NR_SUIT500 + BO_TOTAL) R(NR_SUIT500 + BO_TOTAL)
	- STUDY_AREA covariate omitted from psi and R to make marginal plots for other covariates
	