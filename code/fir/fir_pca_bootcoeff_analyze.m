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
savedir_base = fullfile(masterdir,'analyses','fir');
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
ncomps=8;
%% load bootstrapped components for specified component design
% save with a generic name so the same script can bootstrap resample
%component_design = 'ThreatNonthreatAllStimuliStratified';
savedir = fullfile(savedir_base,['cpc_timecourse_fin',num2str(fin),'st',num2str(st)],component_design,'pncvs22qcoeff');

% reference coefficients; attempt to match order and sign of bootstrapped
% coefficients to these, given arbitrary reordering
OverallCoeffs = load(fullfile(savedir,['FIRGroup',component_design,'_CPCAComponents.mat']));
ref_coeff = OverallCoeffs.nodeDataAll; 
dfun = @(x,y) corr(x,y).^2; % use squared correlation because coefficient sign can arbitrarily flip

d = dir([savedir,'/bootcoeff_AllSubjects/BootstrapPCA_Rep*.mat']);
nreps = length(d);
coeff_boot = nan(nparc_all,ncomps,nreps); % hold all bootstrapped components
coeff_sim = nan(ncomps,nreps); % hold raw similarity between ordinal components
match1 = [4 5 6];
match2 = [6 7 8];
for rep = 1:nreps
    fprintf('Rep %d\n',rep);
    X = load(fullfile(d(rep).folder, d(rep).name)); % load pca results from bootstrapped sample
    % first reorder coefficients based on shared variance, using a limited
    % subset that tend to get flipped
    if CONTAINS(component_design,'IPR')
        try
            [~,shuff_idx] = COLUMN_MATCH(X.coeff(:,match1),ref_coeff(:,match1),1,dfun);
            X.coeff(:,match1) = X.coeff(:,match1(shuff_idx));
        catch
            disp('error: couldnt match components')
        end
        try
            [~,shuff_idx] = COLUMN_MATCH(X.coeff(:,match2),ref_coeff(:,match2),1,dfun);
            X.coeff(:,match2) = X.coeff(:,match2(shuff_idx));
        catch
            disp('error: couldnt match components')
        end
    else
        [~,shuff_idx] = COLUMN_MATCH(X.coeff(:,match1),ref_coeff(:,match1),1,dfun);
        X.coeff(:,match1) = X.coeff(:,match1(shuff_idx));
        [~,shuff_idx] = COLUMN_MATCH(X.coeff(:,match2),ref_coeff(:,match2),1,dfun);
        X.coeff(:,match2) = X.coeff(:,match2(shuff_idx));
    end
    
    coeff_sim(:,rep) = diag(corr(X.coeff(:,1:ncomps),ref_coeff)); % similarity in the current ordering     
    reflect = sign(coeff_sim(:,rep)); % .* abs(coeff_sim(:,rep))>0.95; % now reflect coefficients
    coeff_aligned = X.coeff(:,1:ncomps) .* repmat(reflect',[nparc_all 1]); % reflect coefficients so they all match reference
    coeff_boot(:,:,rep) = coeff_aligned;
end
% compute p(boot weight > 0) and p(boot weight <0)
p = 2*min(cat(3,mean(coeff_boot > 0,3),mean(coeff_boot < 0,3)),[],3); 
p = p*numel(p); 
fprintf('# of components matched at r<0.75: %d\n',sum(sum(abs(coeff_sim)<0.75,1),2));

%% plot similarity distribution
f=figure;
histogram(reshape(abs(coeff_sim),1,[]),'FaceColor',[0.2 0.6 0.3],'FaceAlpha',0.6);
xlabel('Component Loading Similarity (Pearson r)'); ylabel('# of Components');
xlim([0.5 1]);
prettifyEJC; 

%%
surfplot_base = fullfile('analyses','fir',['cpc_timecourse_fin',num2str(fin),'st',num2str(st)],[component_design],'pncvs22qcoeff');
nodeDataAll = OverallCoeffs.nodeDataAll .* (p < 0.05);
nodeData = nodeDataAll(cort_indices,:);
plotTitles = OverallCoeffs.plotTitles;
clim=0.10;
fname_surfplot = ['FIRGroup',component_design,'_CPCAComponentsBootstrappedThreshold.mat'];
save(fullfile(savedir,fname_surfplot),'nodeData','nodeDataAll','plotTitles','clim');

%{
system(['source ~/.bash_profile; source activate pyforge ; python code/visualize/brainvis_fir.py ',name_root,' ',num2str(atlasScale),' ',fullfile(surfplot_base,fname_surfplot)]);            
[img,fname_save] = SURFPLOT_ADD_SUBCORTEX(fullfile(savedir,fname_surfplot),atlasNameSubcortex,atlasScaleSubcortex,sc_indices_combined,clim);
imwrite(img,fname_save);
[dn,fn,~] = fileparts(fname_save);
MPL = load('data/colors/mpl_cmaps.mat');
f=figure; h=colorbar('Location','southoutside');colormap(MPL.custom_ejc1); caxis([-clim clim]);
set(h,'Ticks',[-clim -clim/2 0 clim/2 clim]);set(h,'TickLabels',[-clim -clim/2 0 clim/2 clim])
set(gca,'FontSize',6);
f=FIGURE_SIZE_CM(f,4,4);
saveas(f,fullfile(dn,[fn,'Colormap.pdf']))
%}