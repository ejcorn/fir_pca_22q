% define every possible FIR design matrix

cd(basedir);
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile('data',['TimeSeriesIndicators',name_root,'.mat']));
load(fullfile('data',['ConcTimeSeries',name_root,'.mat']));
masterdir = fullfile('results',name_root);
savedir = fullfile(masterdir,'analyses','fir','design_matrices');
mkdir(savedir);

%% sort out scanids 
[q22mask,demoMatch.is22q] = PROCESS_SCANIDS(demoMatch,subjInd_scanID);

%% set parameters - length of FIR   

TR = 3; nTR = allScanTRs(1);

%% define task response components on all stimuli
allstim = load(fullfile('data/task/idemo/stimulus/all_3col.txt'));
regressor = GET_FIR_REGRESSOR(allstim(:,1),fin,TR,nTR,st);

X = zeros(size(concTS,1),fin*nobs); % make separate regressor for each subject
for N = 1:nobs
    st_x = 1+nTR*(N-1); nd_x = nTR*N;    
    st_y = 1+fin*(N-1); nd_y = fin*N;
    X(st_x:nd_x,st_y:nd_y) = regressor;
end
X = [ones(size(X,1),1) X];
save(fullfile(savedir,['AllStimuli_FIRDesignMatrix_fin',num2str(fin),'st',num2str(st),'.mat']),'X');

%% load correct-incorrect responses
stim_types = {'all','threat','nonthreat'};
stim_types = {'threat','nonthreat'};
for stim_type = stim_types
    stim_type = char(stim_type);
    disp(stim_type);
    CR_stim = load(fullfile(masterdir,['analyses/behavior/idemo/CorrectResponseTiming_',char(stim_type),'.mat']));
    IR_stim = load(fullfile(masterdir,['analyses/behavior/idemo/IncorrectTiming_',char(stim_type),'.mat']));
    NR_stim = load(fullfile(masterdir,['analyses/behavior/idemo/NRTiming_',char(stim_type),'.mat']));
    INR_stim = load(fullfile(masterdir,['analyses/behavior/idemo/IncorrectNRTiming_',char(stim_type),'.mat']));
    regressors = {CR_stim,INR_stim,IR_stim,NR_stim};
    regressor_names = {'correct','incorrect_nr','incorrect','nr'};

    %% tile FIR basis and perform regression

    for i_reg = 1:length(regressor_names)
        disp(regressor_names{i_reg})
        X = zeros(size(concTS,1),fin*nobs); % make separate regressor for each subject
        % indicate whether subject is missing response data or excluded b/c of high NR rate. this is
        % different than having none of a particular response
        
        % subjects missing response data or excluded -> no regressors + exclude BOLD
        % subjects with none of response -> only exclude those regressors
        % set entire regressor column to NaN if you have either of these cases
        % then BOLD data will be excluded based on which rows have nothing
        % in the regressor columns. this SHOULD only correspond to entire
        % subjects who are missing all their response data
        
        SubjectIsMissingOrExcluded = false(nobs,1);
        SubjectNumResponses = nan(nobs,1); % count number of responses per subject
        columnLabelsSubject = nan(1,nobs*fin); % label columns by subject
        columnLabelsTime = nan(1,nobs*fin); % label columns by time point modeled
        
        for N = 1:nobs
            st_x = 1+nTR*(N-1); nd_x = nTR*N;    
            st_y = 1+fin*(N-1); nd_y = fin*N;            
            subj_regressor = regressors{i_reg}.indicator.(['s',subjInd_scanID{1}{N}]); % subject-specific correct-incorrect timing
            columnLabelsTime(st_y:nd_y) = 1:fin; % index time in regressor to extract later
            columnLabelsSubject(st_y:nd_y) = N; % index subject in regressor to extract later
            if strcmp(subj_regressor,'missing') || strcmp(subj_regressor,'excluded')
                SubjectIsMissingOrExcluded(N) = true; 
                % set entire regressor column to NaN if subj missing data
                X(:,st_y:nd_y) = NaN;
            else
                if ~any(isnan(subj_regressor))
                    X(st_x:nd_x,st_y:nd_y) = GET_FIR_REGRESSOR(subj_regressor,fin,TR,nTR,st);
                    SubjectNumResponses(N) = length(subj_regressor);
                else
                    % set entire regressor column to NaN
                    % if subject has none of a particular response
                    X(:,st_y:nd_y) = NaN;
                    disp([subjInd_scanID{1}{N},'-',num2str(N)]);
                end
            end
        end
        
        X = [ones(size(X,1),1) X];    % add intercept
        columnLabelsTime = [0 columnLabelsTime];
        columnLabelsSubject = [0 columnLabelsSubject];
        save(fullfile(savedir,[stim_type,regressor_names{i_reg},'_FIRDesignMatrix_fin',num2str(fin),'st',num2str(st),'.mat']),...
            'X','SubjectIsMissingOrExcluded','SubjectNumResponses','columnLabelsSubject','columnLabelsTime');
    end
end

%% concatenate correct-incorrect, threat-nonthreat into one big predictor matrix

stim_types = {'threat','nonthreat'};
regs = {'correct','incorrect','nr'};
all_reg_data = cell(length(stim_types),length(regs));
for stim_type = 1:length(stim_types)
    for reg = 1:length(regs)
        all_reg_data{stim_type,reg} = load(fullfile(savedir,[char(stim_types(stim_type)),char(regs(reg)),'_FIRDesignMatrix_fin',num2str(fin),'st',num2str(st),'.mat']));
        all_reg_data{stim_type,reg}.stim = char(stim_types(stim_type));
        all_reg_data{stim_type,reg}.reg = char(regs(reg));
    end
end

all_reg_data = reshape(all_reg_data',1,[]);
all_reg_matrix = cellfun(@(X) X.X(:,2:end),all_reg_data,'UniformOutput',false);
X = horzcat(all_reg_matrix{:});

%% label columns
columnLabelsStim = cellfun(@(X) cellstr(repmat(X.stim,[size(X.X,2)-1 1]))',all_reg_data,'UniformOutput',false);
columnLabelsStim = [{'intercept'},columnLabelsStim{:}];
columnLabelsResp = cellfun(@(X) cellstr(repmat(X.reg,[size(X.X,2)-1 1]))',all_reg_data,'UniformOutput',false);
columnLabelsResp = [{'intercept'},columnLabelsResp{:}];
columnLabelsSubject = cellfun(@(X) X.columnLabelsSubject(2:end),all_reg_data,'UniformOutput',false);
columnLabelsSubject = [0,columnLabelsSubject{:}];
columnLabelsTime = cellfun(@(X) X.columnLabelsTime(2:end),all_reg_data,'UniformOutput',false);
columnLabelsTime = [0,columnLabelsTime{:}];
%% check that every time point is modeled or no time points are modeled if subj is missing response data
% now iterate through subjects and make sure that each subject who is not
% missing response data has every time point modeled, which should be
SubjectMissingResponse = cellfun(@(X) X.SubjectIsMissingOrExcluded,all_reg_data,'UniformOutput',false);
SubjectMissingResponse = any(horzcat(SubjectMissingResponse{:}),2); 
all_stim_regressor = GET_FIR_REGRESSOR(allstim(:,1),fin,TR,nTR,st);
SubjectRegressorCount = zeros(nobs,1);
for N = 1:nobs
    st_x = 1+nTR*(N-1); nd_x = nTR*N;  
    SubjectRegressorCount(N) = nansum(nansum(X(st_x:nd_x,:),1),2);
end

if all(SubjectRegressorCount(~SubjectMissingResponse) == sum(sum(all_stim_regressor,1),2))
    disp('Every subject with response data has all stimuli captured by regressors. You are the best eli wow you are so good at matlab it is crazy')
else
    disp('ERROR: some subjects do not have all stimuli captured by regressor');
end
%% test whether matrix is rank deficient and correct

[dup_mask_all] = FIR_REGRESSOR_RM_DUPS(X);

X(:,dup_mask_all) = NaN; % duplicate columns mean solitary stimuli with staggered overlap that cannot be uniquely modeled 

%% add intercept and save

X = [ones(size(X,1),1) X];
save(fullfile(savedir,['ThreatNonthreatAllStimuliStratified_FIRDesignMatrix_fin',num2str(fin),'st',num2str(st),'.mat']),...
    'X','columnLabelsStim','columnLabelsSubject','columnLabelsTime','columnLabelsResp');