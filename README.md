# fir_pca_22q
 
Code to reproduce all analysis in Cornblath et al. 2021 (["Altered functional brain dynamics in chromosome 22q11.2 deletion syndrome during facial affect processing"](https://www.biorxiv.org/content/10.1101/2020.12.17.423342v1)).

## General overview

I have included a basic example of how to apply the CPCA method to fMRI data across subjects in `example.m` in the main folder. This may provide an easier starting point. There is also a publicly available tool box for CPCA.

The script `main_fir.sh` will carry out the sequential execution of all necessary scripts by submitting serially dependent jobs to the Sun Grid Engine job scheduler. 

An HPC is only really necessary for the bootstrapping step to determine significance of the spatial principal component loading. A separate file, `local_main.m` will execute all of the other scripts on your local machine. 

Basic workflow for generating main text figures is as follows:

**Figure 2:**

Scripts:
  - `code/process/ProcessData22q.m`: process BOLD data, demographic data, store as variables that are easy to load in future scripts (`concTS`,`demoMatch`,`subjInd`)
  - `code/behavior/extract_TRlocked_responses_v2.R`: process task response data, i.e. indicators of timing of correct/incorrect responses to each stimulus type
  - `code/fir/fir_design.m`: construct finite impulse response (FIR) design matrices, which are essentially indicator variables for every single task event (combinations of correct/incorrect threat/non-threat stimuli) separately for each subject
  - `code/fir/fir_pca_regression`: perform initial FIR regression, perform PCA on predicted values, then regress component scores (i.e. loading of each PC at each imaging acquisition) back onto FIR design matrix. Stores those regression betas in a matrix.
  - `code/fir/fir_pca_lme_stim.R`: perform multilevel modeling of FIR betas
  - `code/fir/plot_fir_pca_lme_stim.R`: plot results of multilevel modeling
  - `code/fir/fir_pca_bootstrap_prep.m`: prepare subsampled BOLD timeseries and subsampled design matrix for bootstrapping
  - `code/fir/boot_fir_pca.py`: perform FIR CPCA procedure on bootstrapped samples, implemented in python using `sklearn`,`scipy`, and `numpy`. This was done to maximize speed of bootstrapping on HPC.
  - `code/fir/fir_pca_bootcoeff_analyze.m`: 
  - `code/fir/plot_fir_coeffs.m`: make surface plots of FIR PCA spatial components. This is done with the `code/visualize/brainvis_fir.py` script that calls the module `code/visualize/ejcbrain`, which is a wrapper for [pysurfer](https://pysurfer.github.io/). These scripts must be run locally, mostly because I could not get [`mayavi`](https://docs.enthought.com/mayavi/mayavi/) to run on the cluster.

Key Functions:
  - `GET_FIR_REGRESSOR`: create FIR regressor for 1 scan from 1 subject given number of time points, TR, vector of stimulus onsets in seconds
  - `FIR_PCA`: perform initial FIR regression, then perform PCA on the predicted values of that regression, given BOLD data and design matrix
  - `CPCA_SCORES_FIR_REGRESS`: regress PC scores back onto FIR design matrix, given PC scores and design matrix

**Figure 3:**  

Scripts:
  - `code/cognition/fir_pca_cognition_lm.R`: regress FIR betas (or summary metrics of them) onto task performance while adjusting for confounds

**Figure 4:**

Spin test was performed using code from Dr. Aaron Alexander-Bloch's github repository [`spin-test`](https://github.com/spin-test/spin-test). Custom functions in the scripts below are wrappers on functions from the `spin-test` repository.

Scripts:
  - `code/t1/make_spins.m`: generates rotations and is already run for fsaverage5. I've included the data in the repository.
  - `code/t1/t1_spin_pcs.m`: loads the rotations and applies them to difference maps from [Sun et al. 2018](https://www.nature.com/articles/s41380-018-0078-5) (thank you Dr. Sun for generously providing these maps!!)
  - `code/t1/plot_t1_spins.R`: performs non-parametric testing to evaluate the alignment between structural difference and functional spatial PC maps.

## Requirements:
  - BOLD fMRI and T1 MRI data processed using [`fmriprep`](https://fmriprep.org/en/stable/) and [`xcpEngine`](https://xcpengine.readthedocs.io/)
  - MATLAB R2017A or later
  - R 3.6.0 or later, with packages listed in `code/misc/packages.R`    
  - Hardware: a computing cluster using Sun Grid Engine job scheduler, ability to request cores with at least 16G of RAM
  
This software was tested on GNU Linux using the Center for Biomedical Image Computing & Analytics CUBIC computing cluster (https://www.med.upenn.edu/cbica/cubic.html).

## Directory structure and path specification

Scripts are organized by their purpose in the code folder. The jobs folder contains shell scripts to allow the scripts in code folder to be submitted to a computing cluster using a Sun Grid Engine job scheduler (qsub). 

`main_fir.sh` is a bash script that will run all analyses in the paper. There are a few paths that need to be specified in this script:
  - `BASEDIR`: the path to the master branch of this repository, e.g. `/Users/Eli/22qProject/`
  - `MATPATH`: the path to the user's MATLAB binary, e.g. /Applications/MATLAB_R2017b.app/bin/matlab
  - `MATPATHFAST`: can use an older version of MATLAB here to save RAM (we used 2014B), or same as `MATPATH`
  - `RPATH`: the path to the user's Rscript function e.g. Rscript
  - `PYPATH`: the path to user's python binary
  - `LOGDIR`: directory for output files

Additionally, data dependencies are as follows (where `$NPARC` is the parcellation scale and `$ID` is a subject ID):
  - `datadir_main`: at the top of processDataBaumSample.m, specifies the path to processed data. Within this directory there should be:
      - `subject_demographics/*.csv` with csv file of demographics, containing scanids and bblids
      - `diffusion_data/volNormSC/Lausannne${NPARC}/${ID}_volNormStreamline_LausanneScale${NPARC}.mat`
      - `diffusion_data/FA/Lausanne${NPARC}/${ID}_FA_LausanneScale${NPARC}.mat`
      - `rest_data/Lausanne${NPARC}/${ID}_Lausanne${NPARC}_ts.1D`
      - `nback_data/Lausanne${NPARC}/${ID}_Lausanne${NPARC}_ts.1D`
  
`ProcessData_22q.m` requires file paths to the emotion identification task fMRI data processed by `xcpEngine`. To demo this code without obtaining the necessary BOLD data, one could replace the variables "concTS" in this script with random numbers or your own BOLD data.

`code/miscfxns/addpaths.m` also requires specification of path to the BCT:
  - [Brain Connectivity Toolbox](https://sites.google.com/site/bctnet/)  
  - [NIFTI toolbox](https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)

## Input specification

You have the option of specifying a few parameters in `main_fir.sh`:
  - `XCP`: name of XCP output folder
  - `ATLAS`: character name of desired atlas. Must move the nifti file defining your volumetric atlas into `data/nifti` and have XCP output folder matching that parcellation
  - `NPARC`: parcellation scale, i.e. number of nodes in parcellation
  - `SCAN`: character, doesn't matter what you put here. deprecated from original purpose.
  - `LAB`: optional extra label
  - `ST`: how many TRs after stimulus onset to include in your FIR design. 0 means start with the TR in which stimulus onsets, e.g. t= 0 to t = TR
  - `FIN`: how many time points after `ST` to include, e.g. if `ST=0` then you model from t=0 to t=TR\*FIN
  - `ZDIM`=1, only ever set this to 1. Ensures that each region's time series is z-scored within each subject.

In my experience, the entire pipeline takes 3-24 hrs to run. I've never officially timed it, just set and forget :)

## Understanding individual scripts in context

Most of the main analysis scripts are meant to load output from upstream scripts. The following variables are passed into almost every MATLAB script:

 - `basedir`: home/working directory for the project see
 - `name_root`: name of specific output folder for given set of parameters 
 - `st`: number of time points after stimulus to begin modeling BOLD response (0 in paper)
 - `fin`: number of time points after `st` to include in model (6 in paper)
 - `component_design`: specification of design matrix (always "ThreatNonthreatAllStimuliStratified" then may have "\_" followed by additional string specifying a null model to use)

You can access these variables by executing the script `code/miscfxns/STARTUP.m`. If you run this script first, then you can run any other MATLAB script in isolation. Similarly `code/miscfxns/startup_local.R` loads these variables into R.

Please contact Eli Cornblath (Eli ~`DOT`~ Cornblath ~`AT`~ pennmedicine.upenn.edu) with any questions regarding this code.
