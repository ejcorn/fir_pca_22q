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
extralab = extractAfter(name_root,num2str(atlasScale));
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

%% define task response components on all stimuli
%stim_type = 'All'; response_type = 'Stimuli';
%component_design = [component_design];
FIR_design = load(fullfile(savedir_base,'design_matrices',[component_design,'_FIRDesignMatrix_fin',num2str(fin),'.mat']));
X = FIR_design.X;
[coeff_all,scores_all,explained_all,~,betaFIR] = FIR_PCA(X,concTS);

savedir = fullfile(savedir_base,'cpc_timecourse',component_design,'pncvs22qcoeff');
mkdir(savedir);
%% fit separate models
subjInds22q = unique(subjInd(~~SubjectIs22q));
subjIndsPNC = unique(subjInd(~SubjectIs22q));
X_22q_y = [1 find(ismember(FIR_design.columnLabelsSubject,subjInds22q))];
X_22q = X(~~SubjectIs22q,X_22q_y); % get FIR design matrix for just 22q
concTS_22q = concTS(~~SubjectIs22q,:); % get time series matrix for 22q
[coeff_22q,scores_22q,explained_22q] = FIR_PCA(X_22q,concTS_22q); % FIR -> PCA on task related variance in 22q only

X_PNC_y = [1 find(ismember(FIR_design.columnLabelsSubject,subjIndsPNC))];
X_PNC = X(~SubjectIs22q,X_PNC_y); % get FIR design matrix for just PNC
concTS_PNC = concTS(~SubjectIs22q,:); % get time series matrix for PNC
[coeff_PNC,scores_PNC,explained_PNC] = FIR_PCA(X_PNC,concTS_PNC); % FIR -> PCA on task related variance only in PNC only

%% compare models
ncomps=6; ncomps_scree = 50;
PC_idx = 1:ncomps;
dfun = @(x,y) corr(x,y).^2; % using R^2 because maps can be flipped
MPL = load('data/colors/mpl_cmaps.mat');
cmap = MPL.Blues;
f=figure;
subplot(1,4,1);
plot(explained_22q(1:ncomps_scree),'r-'); hold on;
plot(explained_PNC(1:ncomps_scree),'b-'); axis square;
legend({'22q','PNC'});
xlabel('# of Components'); ylabel('Total R^2');
prettifyEJC;

subplot(1,4,2);
[coeffSimilarity,shuffleIdx] =COLUMN_MATCH(coeff_22q(:,PC_idx),coeff_PNC(:,PC_idx),2,dfun);
imagesc(coeffSimilarity); colormap(cmap);
caxis([0 1]); h=colorbar;title(h,'r^2');
xlabel('PNC'); ylabel('22q'); axis square;
set(gca,'FontSize',8);

subplot(1,4,3); % retain indices to rearrange PCs to common order
[coeffSimilarity,shuffleIdx_22qtoAll] =COLUMN_MATCH(coeff_22q(:,PC_idx),coeff_all(:,PC_idx),1,dfun);
imagesc(coeffSimilarity); colormap(cmap);
caxis([0 1]); h=colorbar;title(h,'r^2');
xlabel('Group'); ylabel('22q'); axis square;
set(gca,'FontSize',8);

subplot(1,4,4); % retain indices to rearrange PCs to common order
[coeffSimilarity,shuffleIdx_PNCtoAll] =COLUMN_MATCH(coeff_PNC(:,PC_idx),coeff_all(:,PC_idx),1,dfun);
imagesc(coeffSimilarity); colormap(cmap);
caxis([0 1]); h=colorbar;title(h,'r^2');
xlabel('Group'); ylabel('PNC'); axis square;
f=FIGURE_SIZE_CM(f,18,9);
set(gca,'FontSize',8);

saveas(f,fullfile(savedir,['PNCvs22q',component_design,'FIRPCCoeffs.pdf']))

%% compute variance explained in 22q data by PNC scores and vice versa
% using bootstrapping to generate confidence intervals

ncomps=6; ncomps_scree = 6;
PC_idx = 1:ncomps;
ts_coeff_groups.PNC.concTS = concTS_PNC;
ts_coeff_groups.PNC.coeff = coeff_PNC;
ts_coeff_groups.q22.concTS = concTS_22q;
ts_coeff_groups.q22.coeff = coeff_22q;
ts_coeff_groups.all.coeff = coeff_all;

nboots = 500; nsamps = size(concTS,1);
explained_cross_group = struct();
for coeff_grp = {'PNC','q22','all'}
    coeff_grp = char(coeff_grp);
    for ts_grp = {'PNC','q22'}
        ts_grp = char(ts_grp);
        nsamps = size(ts_coeff_groups.(ts_grp).concTS,1); % get size of time series for bootstrapping
        explained_cross_group.(coeff_grp).(ts_grp) = nan(nboots,nparc_all); % hold bootstrapped explained variance
        for boot = 1:nboots % bootstrap variance explained by group coefficient maps in this sample
            fprintf('%s Coeffs, %s TS, Bootstrap %d\n',coeff_grp,ts_grp,boot)
            bootsamp = randi(nsamps,[nsamps 1]); % indices of bootstrapped sample
            [~,explained_cross_group.(coeff_grp).(ts_grp)(boot,:)] = ...
                PROJECT_DATA_COEFF(ts_coeff_groups.(ts_grp).concTS(bootsamp,:),ts_coeff_groups.(coeff_grp).coeff);
        end
    end
end
%% save csv files for plots
savedir_exp = fullfile(savedir,'bootexplained'); mkdir(savedir_exp);
for coeff_grp = string({'PNC','q22','all'})
    coeff_grp = char(coeff_grp);
    for ts_grp = string({'PNC','q22'})
        ts_grp = char(ts_grp);
        % save first ncomps components, times 2 in case the later
        % components matter a bit for specificity
        tbl = array2table(explained_cross_group.(coeff_grp).(ts_grp)(:,1:(ncomps*2)),'VariableNames',cellfun(@(x) sprintf('PC%d',x),num2cell(1:(ncomps*2)),'UniformOutput',false));
        writetable(tbl,fullfile(savedir_exp,[char(coeff_grp),'Coeff_',char(ts_grp),'TS_Bootstrap.csv']),'Delimiter',',');
    end
end
