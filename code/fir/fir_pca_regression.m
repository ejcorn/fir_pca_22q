% for various types of CPCA components (obtained by PCA of different
% stimulus-linked subsets of the total task related variance)
% measure how those components are expressed over time in PNC and 22q
% as a function of stimulus type and response type (correct incorrect)

cd(basedir); addpath(genpath('code')); addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile('data',['TimeSeriesIndicators',name_root,'.mat']));
load(fullfile('data',['ConcTimeSeries',name_root,'.mat']));
masterdir = fullfile('results',name_root);
savedir_base = fullfile(masterdir,'analyses','fir');
mkdir(savedir_base);
concTS = THRESHOLD(concTS,zdim);

%% get indices of subcortical structures and load subcort BOLD from harvard oxford
atlasNameSubcortex = 'HarvardOxford'; atlasScaleSubcortex = 112;
[concTS,nparc_all,cort_indices,subcort_indices] = CONCAT_CORT_SUBCORT_BOLD(concTS,atlasNameSubcortex,atlasScaleSubcortex);

%% sort out scanids 

[q22mask,demoMatch.is22q] = PROCESS_SCANIDS(demoMatch,subjInd_scanID);

%% set parameters - length of FIR   

fin=6; st = 0;
TR = 3; nTR = allScanTRs(1);
ncomps = 8; % number of components to analyze
resp_thresh = 2; % set minimum number of responses needed to be included in model
%% specify what part of the task related variance you want to get your PCA loadings from

[component_design_load,null_spec] = NULL_SPEC(component_design); % strip away null specification
FIR_Design_Reg1 = load(fullfile(savedir_base,'design_matrices',[component_design_load,'_FIRDesignMatrix_fin',num2str(fin),'st',num2str(st),'.mat']));

%% Implement null models here, if applicable

if strcmp(null_spec,'IPR')
    %f=figure; subplot(1,2,1); imagesc(concTS(1:204,:));
    rng(0); % just doing one rep so set seed to the same value here and in fir_pca_bootstrap_prep
    concTS = IPR_BOLD_NULL(concTS,subjInd);
    %subplot(1,2,2); imagesc(concTS(1:204,:));
elseif strcmp(null_spec,'UPR')
    %f=figure; subplot(1,2,1); imagesc(concTS(1:204,:));
    concTS = UPR_BOLD_NULL(concTS,subjInd);
    %subplot(1,2,2); imagesc(concTS(1:204,:));
end
%%
%f= PLOT_DESIGN_MATRIX(FIR_Design_Reg1,2,subjInd);
%% carry out three-step FIR regression, PCA, FIR regression
 % if using design matrix for initial PCA that contains all stimuli, use these component weights
[coeff,scores,explained] = FIR_PCA(FIR_Design_Reg1.X,concTS); % use scores from input components
% regress those scores back on the full design matrix used to extract
% task-related variance
betaFIR_PCA_all = CPCA_SCORES_FIR_REGRESS(FIR_Design_Reg1.X,scores(:,1:ncomps)); 

% make a separate folder for analysis of each PCA solution
savedir = fullfile(masterdir,'analyses','fir',['cpc_timecourse_fin',num2str(fin),'st',num2str(st)],component_design);
mkdir(savedir);
save(fullfile(savedir,'GroupCPCAComponentsExplained.mat'),'explained'); % save explained to make scree plot in R
mkdir(fullfile(savedir,'pncvs22qcoeff'));
writetable(cell2table({['PCA performed on task related variance from ',component_design,', those scores regressed back onto various responses.']}),fullfile(savedir,'description.txt'))
%% load correct-incorrect responses
% loop through stimulus-response combinations and 
% extract betas from regression matrix into a more workable format
% for linear mixed effects modeling

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
