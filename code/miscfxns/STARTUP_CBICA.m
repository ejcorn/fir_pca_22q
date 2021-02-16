% start every script for testing on CBICA, specific to a certain name root
close all; clc
basedir = '/cbica/home/cornblae/ecornblath/fir_pca_22q';
name_root = 'CPCA_IDSchaefer200Z1xcp_6p_noFilter'; 
cd(basedir); addpath(genpath('code'));