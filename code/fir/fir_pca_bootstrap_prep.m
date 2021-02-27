% compare FIR components between controls, 22q, and group
% to test whether your temporal findings primarily are driven
% by spatial differences

cd(basedir);
run('code/miscfxns/addpaths.m');
cd(basedir);
addpath(genpath('code'));
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile('data',['TimeSeriesIndicators',name_root,'.mat']));
load(fullfile('data',['ConcTimeSeries',name_root,'.mat']));
masterdir = fullfile('results',name_root);
savedir_base = fullfile(masterdir,'analyses','fir');
mkdir(savedir_base);
concTS = THRESHOLD(concTS,zdim);

%% get indices of subcortical structures and load subcort BOLD from harvard oxford
atlasNameSubcortex = 'HarvardOxford'; atlasScaleSubcortex = 112;
[concTS,nparc_all,cort_indices,sc_indices_combined] = CONCAT_CORT_SUBCORT_BOLD(concTS,atlasNameSubcortex,atlasScaleSubcortex,atlasScale,name_root);

%% sort out scanids 
[q22mask,demoMatch.is22q] = PROCESS_SCANIDS(demoMatch,subjInd_scanID);

%% set parameters - length of FIR   

TR = 3; nTR = allScanTRs(1);
ncomps = 10; % save more components to give some wiggle room for matching components in case amount of variance 
% explained is variable, e.g. PC6 of one sample is PC8 of another sample because PC6-8 all explain similar amount of variance

%% Load design matrix and implement null models here, if applicable

[component_design_load,null_spec] = NULL_SPEC(component_design); % strip away null specification
FIR_Design_Reg1 = load(fullfile(savedir_base,'design_matrices',[component_design_load,'_FIRDesignMatrix_fin',num2str(fin),'st',num2str(st),'.mat']));
[concTS,FIR_Design_Reg1] = NULL_IMPLEMENT(null_spec,concTS,subjInd,FIR_Design_Reg1,name_root); % see functions for descriptions of null models

%% split regressors into 22q half and PNC half
% save with a generic name so the same script can bootstrap resample

FIR_design = load(fullfile(savedir_base,'design_matrices',[component_design_load,'_FIRDesignMatrix_fin',num2str(fin),'st',num2str(st),'.mat']));

savedir = fullfile(savedir_base,['cpc_timecourse_fin',num2str(fin),'st',num2str(st)],[component_design],'pncvs22qcoeff');
mkdir(savedir);

X_all = FIR_design.X; % make it so every time you save variables named concTS and X, even if whole sample or groups
concTS_all = concTS;
X = X_all;
save(fullfile(savedir,['AllSubjectsBOLD',component_design,'Regressor.mat']),'concTS','X','ncomps');

subjInds22q = unique(subjInd(~~SubjectIs22q));
subjIndsPNC = unique(subjInd(~SubjectIs22q));
X_22q_y = [1 find(ismember(FIR_design.columnLabelsSubject,subjInds22q))];
X = X_all(~~SubjectIs22q,X_22q_y); % get FIR design matrix for just 22q
concTS = concTS_all(~~SubjectIs22q,:); % get time series matrix for 22q
save(fullfile(savedir,['22qBOLD',component_design,'Regressor.mat']),'concTS','X','ncomps');

X_PNC_y = [1 find(ismember(FIR_design.columnLabelsSubject,subjIndsPNC))];
X = X_all(~SubjectIs22q,X_PNC_y); % get FIR design matrix for just PNC
concTS = concTS_all(~SubjectIs22q,:); % get time series matrix for PNC

save(fullfile(savedir,['PNCBOLD',component_design,'Regressor.mat']),'concTS','X','ncomps');
