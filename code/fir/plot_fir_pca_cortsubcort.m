% plot FIR components generated by
% fir_pca_correctincorrect_cortsubcort.m


addpaths;
cd(basedir);
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile('data',['TimeSeriesIndicators',name_root,'.mat']));
load(fullfile('data',['ConcTimeSeries',name_root,'.mat']));
masterdir = fullfile('results',name_root);
savedir = fullfile(masterdir,'analyses','fir','subject_fir_correct_incorrect_pca');
mkdir(savedir);
concTS = THRESHOLD(concTS,zdim);
%% get indices of subcortical structures and load subcort BOLD from brainnetome
atlasNameSubcortex = 'Laus'; atlasScaleSubcortex = 125;
[nifti,sc_indices] = RETURN_NII_SUBCORT(atlasNameSubcortex,atlasScaleSubcortex);
extralab = extractAfter(name_root,num2str(atlasScale));
% load time series data from
concTS_subcort = load(fullfile('data',sprintf(['ConcTimeSeriesScanID%s%d',extralab,'.mat'],atlasNameSubcortex,atlasScaleSubcortex)));
%% concatenate brainnetome subcortex with schaefer subcortex
concTS = [concTS concTS_subcort.concTS(:,sc_indices)];
cort_indices = 1:nparc; % indices for schaefer cortical parcels
nparc_all = size(concTS,2);
sc_indices_combined = (nparc+1):nparc_all; % location of subcortical nodes in concatenated matrix

%% set parameters - length of FIR   

fin=6;
TR = 3; nTR = allScanTRs(1);
%% load correct-incorrect responses
stim_types = {'threat','nonthreat'};
regressor_names = {'correct','incorrect_nr'};

fname_surfplot = fullfile('analyses','fir','subject_fir_correct_incorrect_pca','FIRAllStimuli_CPCAComponents.mat');
system(['source ~/.bash_profile; source activate pyforge ; python code/assesscluster/brainvis_fir.py ',name_root,' ',num2str(atlasScale),' ',fname_surfplot]);            

    %}
    %
%% add subcortex to those plots
clim=0.1;
MPL = load('data/colors/mpl_cmaps.mat');
cmap = MPL.cividis;
FIGURE_DISPLAY('off');        
% load data for plot type produced by
fname_in = fullfile(masterdir,fname_surfplot);
plot_data = load(fname_in);
n_comps = size(plot_data.nodeDataAll,2);
%
for comp = 1:n_comps
    fname_subcort = sprintf('FIR%d%s%d',comp,atlasNameSubcortex,atlasScaleSubcortex);

    nodeData = plot_data.nodeDataAll(:,comp); % get node data for specific plot
    [f_ant,f_lat] = plot_subcortvol(nodeData(sc_indices_combined),sc_indices,sc_indices,nifti,cmap,...
        -clim,clim);
    set(0, 'CurrentFigure', f_lat); %colorbar; % give colorbar to lateral
    f_ant = FIGURE_SIZE_CM(f_ant,4,4);
    f_lat = FIGURE_SIZE_CM(f_lat,4,4);
    saveas(f_ant,fullfile(savedir,[fname_subcort,'ant.png']));
    saveas(f_lat,fullfile(savedir,[fname_subcort,'lat.png']));
    CData1 = imread(fullfile(savedir,[fname_subcort,'ant.png']));
    CData2 = imread(fullfile(savedir,[fname_subcort,'lat.png']));
    CDataCombined = HORZ_WHITEPAD([CData1 CData2],1.25*2*size(CData1,2)); % white pad           
    imwrite(CDataCombined,fullfile(savedir,[fname_subcort,'.png']));
    delete(fullfile(savedir,[fname_subcort,'ant.png']));
    delete(fullfile(savedir,[fname_subcort,'lat.png']));
end
%
% join cortex plots with subcortex plots

[~,fname_surfplot,~] = fileparts(fname_surfplot);
if exist(fullfile(savedir,[fname_surfplot,'.png']),'file')
    cortexPlot = imread(fullfile(savedir,[fname_surfplot,'.png']));
    subcortexFiles = cellfun(@(comp) fullfile(savedir,sprintf('FIR%d%s%d.png',comp,atlasNameSubcortex,atlasScaleSubcortex)),num2cell(1:n_comps),'UniformOutput',false);
    subcortexPlots = cellfun(@(f) imread(f),subcortexFiles,'UniformOutput',false);
    subcortexPlots = cat(2,subcortexPlots{:}); % horizontally concatenate all plots into one
    subcortexPlots = imresize(subcortexPlots,size(cortexPlot,2)/size(subcortexPlots,2)); % make same width as surface plots     
    imwrite([cortexPlot;subcortexPlots],fullfile(savedir,[fname_surfplot,'Subcortex.png']));        

    for f = subcortexFiles
         delete(char(f));
    end

end

FIGURE_DISPLAY('on');