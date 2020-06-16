% compare FIR components between controls, 22q, and group
% to test whether your temporal findings primarily are driven
% by spatial differences

cd(basedir);
run('code/miscfxns/addpaths.m');
cd(basedir);
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile('data',['TimeSeriesIndicators',name_root,'.mat']));
load(fullfile('data',['ConcTimeSeries',name_root,'.mat']));
masterdir = fullfile('results',name_root);
savedir_base = fullfile(masterdir,'analyses','fir','subject_fir_correct_incorrect_pca');
mkdir(savedir_base);
concTS = THRESHOLD(concTS,zdim);
%% get indices of subcortical structures and load subcort BOLD from brainnetome
atlasNameSubcortex = 'HarvardOxford'; atlasScaleSubcortex = 112;
[nifti,sc_indices] = RETURN_NII_SUBCORT(atlasNameSubcortex,atlasScaleSubcortex);
extralab = extractAfter_(name_root,num2str(atlasScale));
% load time series data from
concTS_subcort = load(fullfile('data',sprintf(['ConcTimeSeriesCPCA_ID%s%d',extralab,'.mat'],atlasNameSubcortex,atlasScaleSubcortex)));
%% concatenate brainnetome subcortex with schaefer subcortex
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
ncomps = 10;
%% split regressors into 22q half and PNC half
% save with a generic name so the same script can bootstrap resample
FIR_design = load(fullfile(savedir_base,'design_matrices',[component_design,'_FIRDesignMatrix_fin',num2str(fin),'.mat']));

savedir = fullfile(savedir_base,'cpc_timecourse',[component_design],'pncvs22qcoeff');
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