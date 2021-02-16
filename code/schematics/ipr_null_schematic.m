% plot IPR null model for figure as example of what null is doing 
% this top part of script copied from code/fir/fir_pca_regression.m
cd(basedir); addpath(genpath('code')); addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile('data',['TimeSeriesIndicators',name_root,'.mat']));
load(fullfile('data',['ConcTimeSeries',name_root,'.mat']));
masterdir = fullfile('results',name_root);
savedir_base = fullfile(masterdir,'analyses','fir','subject_fir_correct_incorrect_pca');
mkdir(savedir_base);
concTS = THRESHOLD(concTS,zdim);
%% get indices of subcortical structures and load subcort BOLD from harvard oxford
atlasNameSubcortex = 'HarvardOxford'; atlasScaleSubcortex = 112;
[nifti,sc_indices] = RETURN_NII_SUBCORT(atlasNameSubcortex,atlasScaleSubcortex);
extralab = extractAfter_(name_root,num2str(atlasScale));
% load time series data from
concTS_subcort = load(fullfile('data',sprintf(['ConcTimeSeriesCPCA_ID%s%d',extralab,'.mat'],atlasNameSubcortex,atlasScaleSubcortex)));
%% concatenate harvard oxford subcortex with schaefer subcortex
concTS = [concTS concTS_subcort.concTS(:,sc_indices)];
cort_indices = 1:nparc; % indices for schaefer cortical parcels
nparc_all = size(concTS,2);
sc_indices_combined = (nparc+1):nparc_all; % location of subcortical nodes in concatenated matrix

%% sort out scanids 
q22mask = ismember(subjInd_scanID{1},cellstr(num2str(demoMatch.scanid(strcmp(demoMatch.study,'22q')))));
demoMatch.is22q = double(strcmp(demoMatch.study, '22q'));

%% set parameters - length of FIR   

fin=6;
TR = 3; nTR = allScanTRs(1);
ncomps = 8; % number of components to analyze
resp_thresh = 2; % set minimum number of responses needed to be included in model
%% specify what part of the task related variance you want to get your PCA loadings from

[component_design_load,null_spec] = NULL_SPEC(component_design); % strip away null specification
FIR_Design_Reg1 = load(fullfile(savedir_base,'design_matrices',[component_design_load,'_FIRDesignMatrix_fin',num2str(fin),'.mat']));

%% Implement null models here, if applicable
rng(0); % just doing one rep so set seed to the same value here and in fir_pca_bootstrap_prep
concTS_null = IPR_BOLD_NULL(concTS,subjInd);
%% plot null models for 1 subject
% next to convolved stimulus
allstim = load(fullfile('data/task/idemo/stimulus/all_3col.txt'));
spmpath = '~/Dropbox/ToolBoxes/spm12/';
addpath(spmpath);
[conv_stim,boxcar] = CONV_STIM_BOXCAR(allstim,5.5,TR,nTR);

subj = 100;
subj_mask = subjInd == subj;
X = concTS(subj_mask,:);
X_null = concTS_null(subj_mask,:);

ol_scl = 26; % scale convolved stimulus overlay
t_tr = linspace(0,nTR,4); % time in trs
% actual data, all regions
f=figure; subplot(1,3,1); imagesc(X'); hold on;
plot(-ol_scl*conv_stim -max(-ol_scl*conv_stim)+nparc_all,'w','LineWidth',1); 
title('BOLD'); ylabel('Region'); xlabel('Time (s)');
ylim([0 nparc_all]); xlim([0 nTR]);
xticks(t_tr); xticklabels(TR*t_tr);
COLORBAR(round(max(abs(X),[],'all'),2,'significant'))
set(gca,'FontSize',8,'TickLength',[0 0]);

% null data, all regions
subplot(1,3,2); imagesc(X_null'); hold on;
plot(-ol_scl*conv_stim -max(-ol_scl*conv_stim)+nparc_all,'w','LineWidth',1); 
title('Phase Randomized BOLD'); ylabel('Region'); xlabel('Time (s)');
ylim([0 nparc_all]); xlim([0 nTR]);
xticks(t_tr); xticklabels(TR*t_tr);
COLORBAR(round(max(abs(X),[],'all'),2,'significant'))
set(gca,'FontSize',8,'TickLength',[0 0]);

% plot left occipital pole (schaefer 200, parcel 5)
parcel = 5;
v1 = X(:,parcel); v1_null = X_null(:,parcel);
subplot(1,3,3); plot(v1,'r'); hold on; plot(concTS_null(1:204,5),'b');
plot(conv_stim+min(v1)-1,'k'); 
legend({['V1, r=',LABELROUND2(corr(conv_stim,v1))],...
    ['V1 Null, r=',LABELROUND2(corr(conv_stim,v1_null))]});
ylim([min(v1)-1,max(v1)+1]); xlim([0 nTR]);
xticks(t_tr); xticklabels(TR*t_tr);
title('Left Occipital Pole'); ylabel('BOLD (z)'); 
xlabel('Time (s)');
prettifyEJC
f=FIGURE_SIZE_CM(f,21,3)
savedir = fullfile(masterdir,'analyses','fir','subject_fir_correct_incorrect_pca','cpc_timecourse',component_design);
saveas(f,fullfile(savedir,'IPR_NULL_schematic.pdf'));
%% confirm null destroys relationship between stimulus and signal
% distribution of correlations between left v1 signal and convolved
% stimulus, for real data and null
parcel = 5;
results = struct();
results.r = zeros(nobs,1);
results.r_null = zeros(nobs,1);
for subj = 1:nobs
    results.r(subj) = corr(concTS(subjInd==subj,parcel),conv_stim);
    results.r_null(subj) = corr(concTS_null(subjInd==subj,parcel),conv_stim);
end
f=figure; histogram(results.r); hold on; histogram(results.r_null);
