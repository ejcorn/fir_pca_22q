% for various types of CPCA components (obtained by PCA of different
% stimulus-linked subsets of the total task related variance)
% measure how those components are expressed over time in PNC and 22q
% as a function of stimulus type and response type (correct incorrect)

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
ncomps = 8; % number of components to analyze
resp_thresh = 2; % set minimum number of responses needed to be included in model
%% specify what part of the task related variance you want to get your PCA loadings from
%component_design = 'allcorrect';
FIR_Design_Reg1 = load(fullfile(savedir_base,'design_matrices',[component_design,'_FIRDesignMatrix_fin',num2str(fin),'.mat']));
 % if using design matrix for initial PCA that contains all stimuli, use these component weights
[coeff,scores,explained] = FIR_PCA(FIR_Design_Reg1.X,concTS); % use scores from input components
% regress those scores back on the full design matrix used to extract
% task-related variance
betaFIR_PCA_all = CPCA_SCORES_FIR_REGRESS(FIR_Design_Reg1.X,scores(:,1:ncomps)); 

% make a separate folder for analysis of each PCA solution
savedir = fullfile(masterdir,'analyses','fir','subject_fir_correct_incorrect_pca','cpc_timecourse',component_design);
mkdir(savedir);
save(fullfile(savedir,'GroupCPCAComponentsExplained.mat'),'explained'); % save explained to make scree plot in R
mkdir(fullfile(savedir,'pncvs22qcoeff'));
writetable(cell2table({['PCA performed on task related variance from ',component_design,', those scores regressed back onto various responses.']}),fullfile(savedir,'description.txt'))
%% load correct-incorrect responses

stim_types = {'threat','nonthreat'};
for stim_type = stim_types
    stim_type = char(stim_type);    
    regressor_names = {'correct','incorrect'};
    for i_reg = 1:length(regressor_names)
        disp([stim_type,'-',regressor_names{i_reg}]);
        stim_response_mask = strcmp(FIR_Design_Reg1.columnLabelsStim,stim_type) & strcmp(FIR_Design_Reg1.columnLabelsResp,regressor_names{i_reg});
        %% visualize estimated stimulus response       
        betaFIR_reshape = zeros(nobs,fin,ncomps); % reshape betas for statistics
        for tp = 1:fin
            % get beta for ith time point and all PCs for each subject, excluding intercept
            betaFIR_reshape(:,tp,:) = betaFIR_PCA_all(FIR_Design_Reg1.columnLabelsTime==tp & stim_response_mask,1:ncomps); % get betas in subjxcomponent matrix for given time point                  
        end                 
        
        %% reflect coefficient axes (and scores) so that time courses are positive for threat correct
        % find the direction (positive or negative) of each score's time course 
        if strcmp(stim_type,'threat') && strcmp(regressor_names{i_reg},'correct') % this is the first iteration and reflects the PCs based upon that
            score_signs = sign(nanmean(nanmean(betaFIR_reshape(q22mask==0,:,:),2),1)); % this is a stick poking back into the 3rd dimension
            % reflect coefficients 
            coeff_reflect = coeff(:,1:ncomps) .* repmat(reshape(squeeze(score_signs),1,[]),[nparc_all 1]);
        end
        % reflect scores - tile that stick (score_signs) in rows and columns
        betaFIR_reshape = betaFIR_reshape .* repmat(score_signs,[nobs fin 1]);
        %% save scores for R to analyze time courses with mixed effects models
        
        save(fullfile(savedir,[stim_type,regressor_names{i_reg},'FIRBetas_CPCScores.mat']),'betaFIR_reshape','subjInd_scanID','explained');
        
    end
end

%% save coefficients to plot later using variable names for brainvis_nodedata.py 
nodeDataAll = coeff_reflect;
explained_tmp = round(explained(1:ncomps),1,'decimals')';
nodeData = nodeDataAll(cort_indices,:); 
plotTitles = cellfun(@(x,y) sprintf('PC%d: %0.1f%%',x,y),GET_PC_IDX(ncomps,XCP_folder),num2cell(explained_tmp),'UniformOutput',false)';
clim=0.10;
fname_surfplot = ['FIRGroup',component_design,'_CPCAComponents.mat'];
save(fullfile(savedir,'pncvs22qcoeff',fname_surfplot),'nodeData','nodeDataAll','plotTitles','clim');