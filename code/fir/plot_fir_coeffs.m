% plot cortex and subcortex of CPC coefficient spatial maps

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
[concTS,nparc_all,cort_indices,sc_indices_combined] = CONCAT_CORT_SUBCORT_BOLD(concTS,atlasNameSubcortex,atlasScaleSubcortex,atlasScale,name_root);

%%
% specify component design, corresponding to which responses are modeled
%component_design = 'ThreatNonthreatAllStimuliStratified';
savedir = fullfile(savedir_base,['cpc_timecourse_fin',num2str(fin),'st',num2str(st)],component_design,'pncvs22qcoeff');
surfplot_base = fullfile('analyses','fir',['cpc_timecourse_fin',num2str(fin),'st',num2str(st)],[component_design],'pncvs22qcoeff');

% plot both thresholded and unthresholded maps
%fnames_surfplot = {['FIRGroup',component_design,'_CPCAComponents.mat'],['FIRGroup',component_design,'_CPCAComponentsBootstrappedThreshold.mat']};
fnames_surfplot = {['FIRGroup',component_design,'_CPCAComponents.mat']};
%fnames_surfplot = {['FIRGroup',component_design,'_CPCAComponentsBootstrappedThreshold.mat']};

for fname_surfplot = fnames_surfplot
	fname_surfplot = char(fname_surfplot);
	load(fullfile(savedir,fname_surfplot),'clim'); % load color limits for that plot

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
end