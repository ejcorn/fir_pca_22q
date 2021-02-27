% start every script for local testing, specific to a certain name root
close all; clc
basedir = '~/Dropbox/Cornblath_Bassett_Projects/BrainStates22q/fir_pca_22q';
name_root = 'CPCA_IDSchaefer200Z1xcp_6p_noFilter'; 
cd(basedir); addpath(genpath('code'));
fin = 6; st = 0; 