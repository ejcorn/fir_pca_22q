%% local testing

%% fir design matrices
clear all; component_design = 'ThreatNonthreatAllStimuliStratified'; STARTUP;
run('code/fir/fir_pca_regression.m');
%%
clear all; component_design = 'ThreatNonthreatAllStimuliStratified'; STARTUP;
run('code/fir/fir_pca_bootcoeff_analyze.m');
%%
clear all; component_design = 'ThreatNonthreatAllStimuliStratified'; STARTUP;
run('code/fir/plot_fir_coeffs.m');

%%
run('code/fir/fir_design.m');
%% run FIR pca all
%component_designs = {'allincorrect','allcorrect','threatcorrect','threatincorrect','nonthreatcorrect','nonthreatincorrect'};
component_designs = {'ThreatNonthreatAllStimuliStratified'};%,'allincorrect','allcorrect'};
for component_design = component_designs
    component_design = char(component_design);
    clearvars -except component_design; STARTUP;
    run('code/fir/fir_pca_correctincorrect_cortsubcort_R.m');
    if strcmp(component_design,'ThreatNonthreatAllStimuliStratified')
        clearvars -except component_design; STARTUP;
        %run('code/fir/fir_pca_compare_coeffs_ci.m');
    end
end
%%
clear all; close all; clc
basedir = '~/Dropbox/Cornblath_Bassett_Projects/BrainStates22q/brain_states_22q';
name_root = 'ScanIDSchaefer200Z1_22q_PreQC_XCP36pdespike_100reps'; 
cd(basedir); addpath(genpath('code'));
numClusters = 6;
run('code/fir/fir_metrics_ci_cortsubcort.m');

clear all; close all; clc
basedir = '~/Dropbox/Cornblath_Bassett_Projects/BrainStates22q/brain_states_22q';
name_root = 'ScanIDSchaefer200Z1_22q_PreQC_XCP36pdespike_100reps'; 
cd(basedir); addpath(genpath('code'));
numClusters = 6;
run('code/fir/plot_fir_metrics_ci_cortsubcort.m');
%% brainnetome
clear all; close all; clc
basedir = '~/Dropbox/Cornblath_Bassett_Projects/BrainStates22q/brain_states_22q';
name_root = 'ScanIDBrainnetome246Z1_22q_100reps'; 
%basedir = '/data/tesla-data/ecornblath/brain_states_22q/'; 
cd(basedir); addpath(genpath('code'));