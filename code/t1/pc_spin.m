addpaths;
cd(basedir);
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile('data',['TimeSeriesIndicators',name_root,'.mat']));
masterdir = fullfile('results',name_root);
savedir_base = fullfile(masterdir,'analyses','t1');
mkdir(savedir_base);

%% load components

component_design = 'ThreatNonthreatAllStimuliStratified';
savedir_func = fullfile(masterdir,'analyses','fir',...
    ['cpc_timecourse_fin',num2str(fin),'st',num2str(st)],component_design,'pncvs22qcoeff');
fname_components = ['FIRGroup',component_design,'_CPCAComponentsBootstrappedThreshold.mat'];
load(fullfile(savedir_func,fname_components),'nodeData','nodeDataAll','plotTitles','clim');

% make final save directory
savedir_comp = fullfile(savedir_base,'pc_spin',component_design);
mkdir(savedir_comp);
%% load rotations of fsaverage5 surface
subject = 'fsaverage5';

% spin a bunch of times and save the rotation
load(fullfile(savedir_base,['SpinTestRotations',subject,'.mat']),'rotateRight','rotateLeft');

nperms = 500;

%% apply rotations to data
% load precomputed spins
testMaps = [nodeData]; % test alignment between PCs and T1 morphometry
% convert centroid difference maps to fsaverage5 surface
[lh_surf_data,rh_surf_data] = SCHAEFER_NODE_TO_SURF_FAST(testMaps,subject);

% apply rotations
lh_surf_data_perm = nan([size(lh_surf_data) nperms]);
rh_surf_data_perm = nan([size(rh_surf_data) nperms]);
testMaps_perm = nan([size(testMaps) nperms]);

for p = 1:nperms
    sprintf('Permutation %d',p)
    lh_surf_data_perm(:,:,p) = lh_surf_data(rotateLeft{p},:); % rotate all data simultaneously with each stored rotation
    rh_surf_data_perm(:,:,p) = rh_surf_data(rotateRight{p},:); % rotate all data simultaneously with each stored rotation
    testMaps_perm(:,:,p) = SCHAEFER_SURF_TO_NODE_FAST(lh_surf_data_perm(:,:,p),rh_surf_data_perm(:,:,p),nparc,subject);
end

%% load reference maps -- T1 structural data from Sun et al.
% load betas
d = dir('data/Sunetal2018_Fig1Maps/*beta*.mgh');
fnames = strcat({d.folder}','/',{d.name}');
mgh_data = cellfun(@(x) load_mgh(x),fnames,'Uni',false);
mgh_data = cellfun(@(x) squeeze(x(:,:,:,2)),mgh_data,'UniformOutput',false);
t1_data = [mgh_data{:}];
[ct_data] = SCHAEFER_SURF_TO_NODE_FAST(t1_data(:,1),t1_data(:,3),nparc,'fsaverage');
[sa_data] = SCHAEFER_SURF_TO_NODE_FAST(t1_data(:,2),t1_data(:,4),nparc,'fsaverage');

save(fullfile(savedir_comp,['SpinTestComponents',component_design,'.mat']),'testMaps','testMaps_perm','ct_data','sa_data');
%% null vs actual variance explained

cor_type='spearman';
refmap_cell = {ct_data,sa_data};
refmap_name = {'Cortical Thickness','Cortical Surface Area'};
simfun = @(t1,bold) corr(t1,abs(bold),'rows','pairwise','type',cor_type).^2;
r2_actual = cell(length(refmap_cell),1);
r2_null = cell(length(refmap_cell),1);
p_spin = cell(length(refmap_cell),1);
p_pearson = cell(length(refmap_cell),1);

% don't compute correlations over 0 areas
%testMaps(testMaps==0) = NaN;
%testMaps_perm(testMaps_perm ==0) = NaN;

for j = 1:length(refmap_cell)
    r2_actual{j} = simfun(refmap_cell{j},testMaps);
    r2_null{j} = nan(nperms,size(testMaps,2));
    for p = 1:nperms
        r2_null{j}(p,:) = simfun(refmap_cell{j},testMaps_perm(:,:,p));
    end
    % p value interpretation: probability that you observe greater R^2 in
    % null
    p_spin{j} = mean(r2_actual{j} < r2_null{j},1);
    [~,p_pearson{j}] = corr(refmap_cell{j},testMaps);
end

%% plot

n_refs = length(refmap_cell);
n_PCs = size(testMaps,2);

msz=20; nsig = 2; l = lines; mcol = l(5,:);
for ref = 1:n_refs
    refMap_i = refmap_cell{ref};    
    f=figure; [sp1,sp2] = subplot_ind2(n_PCs);
    for PC = 1:n_PCs        
        %subplot(n_refs,n_PCs,(ref-1)*n_PCs+PC);
        subplot(sp2,sp1,PC);
        scatter(refMap_i,abs(testMaps(:,PC)),msz,'filled','MarkerFaceColor',mcol,'MarkerFaceAlpha',0.5); hold on;
        lsline; axis square;
        xlabel(refmap_name{ref}); ylabel(sprintf('PC%d',PC)); 
        %title({['R^2 = ',LABELROUND2((r2_actual{ref}(PC)),nsig)],[' p = ',LABELROUND2(p_pearson{ref}(PC),nsig)],...
            %['p_{spin} = ',LABELROUND2(p_spin{ref}(PC),nsig)]});
        title(['R^2 = ',LABELROUND2((r2_actual{ref}(PC)),nsig),', p = ',LABELROUND2(p_pearson{ref}(PC),nsig),...
            ', p_{spin} = ',LABELROUND2(p_spin{ref}(PC),nsig)]);
        prettifyEJC;
    end
    f=FIGURE_SIZE_CM(f,20,9);
    saveas(f,fullfile(savedir_comp,['PCsVs',refmap_name{ref},cor_type,'.pdf']));
end

%%
sqrt(vertcat(r2_actual{:}))
vertcat(p_spin{:})
f=figure;
plot(refmap_cell{2},testMaps(:,5),'b.');
