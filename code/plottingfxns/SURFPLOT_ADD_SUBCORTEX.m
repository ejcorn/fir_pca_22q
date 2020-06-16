function [img,fname_save] = SURFPLOT_ADD_SUBCORTEX(fname_surfplot,atlasNameSubcortex,atlasScaleSubcortex,sc_indices_combined,clim,cmap)
    % INPUTS:
    % fname_surfplot: name of .mat file containing a variable 'nodeDataAll'
    %   -nodeDataAll is an NxK matrix  that has values associated with
    %   parcels
    %   a file with same name as fname_surfplot but .png extension containing
    %   cortical surface plots should exist
    % atlasNameSubcortex: name of atlas used for subcortical parcels in
    % nodeDataAll
    % atlasScaleSubcortex: scale of atlas used for subcortical parcels
    % sc_indices_combined: location of subcortical nodes in nodeData
    % clim: symmetric color limit
    % cmap: colormap
    %
    % OUTPUTS:
    % img: bitmap image merging surface plot and anterior + lateral
    % subcortex volume meshes
    % fname_save: fname_surfplot file name + Subcortex + .png
    % use imwrite(img,fname_save) outside this function to save file
    
    [nifti,sc_indices] = RETURN_NII_SUBCORT(atlasNameSubcortex,atlasScaleSubcortex); % load nifti file for subcortical atlas. requires files to exist in relative directory    
    MPL = load('data/colors/mpl_cmaps.mat');
    if ~exist('cmap','var')
        cmap = MPL.custom_ejc1(:,1:3);
    end
    [folder_surfplot,fname_surfplot,~] = fileparts(fname_surfplot);

    FIGURE_DISPLAY('off');        
    % load data for plot produced by
    plot_data = load(fullfile(folder_surfplot,[fname_surfplot,'.mat']));
    K = size(plot_data.nodeDataAll,2);
    %
    for j = 1:K
        fname_subcort = sprintf('Subcortex%d%s%d',j,atlasNameSubcortex,atlasScaleSubcortex);

        nodeData = plot_data.nodeDataAll(:,j); % get node data for specific plot
        [f_ant,f_lat] = plot_subcortvol(nodeData(sc_indices_combined),sc_indices,sc_indices,nifti,cmap,...
            -clim,clim);
        set(0, 'CurrentFigure', f_lat); %colorbar; % give colorbar to lateral
        f_ant = FIGURE_SIZE_CM(f_ant,4,4);
        f_lat = FIGURE_SIZE_CM(f_lat,4,4);
        saveas(f_ant,fullfile(folder_surfplot,[fname_subcort,'ant.png']));
        saveas(f_lat,fullfile(folder_surfplot,[fname_subcort,'lat.png']));
        CData1 = imread(fullfile(folder_surfplot,[fname_subcort,'ant.png']));
        CData2 = imread(fullfile(folder_surfplot,[fname_subcort,'lat.png']));
        CDataCombined = HORZ_WHITEPAD([CData1 CData2],1.25*2*size(CData1,2)); % white pad           
        imwrite(CDataCombined,fullfile(folder_surfplot,[fname_subcort,'.png']));
        delete(fullfile(folder_surfplot,[fname_subcort,'ant.png']));
        delete(fullfile(folder_surfplot,[fname_subcort,'lat.png']));
    end
    %
    % join cortex plots with subcortex plots

    if exist(fullfile(folder_surfplot,[fname_surfplot,'.png']),'file')
        cortexPlot = imread(fullfile(folder_surfplot,[fname_surfplot,'.png']));
        subcortexFiles = cellfun(@(comp) fullfile(folder_surfplot,sprintf('Subcortex%d%s%d.png',comp,atlasNameSubcortex,atlasScaleSubcortex)),num2cell(1:K),'UniformOutput',false);
        subcortexPlots = cellfun(@(f) imread(f),subcortexFiles,'UniformOutput',false);
        subcortexPlots = cat(2,subcortexPlots{:}); % horizontally concatenate all plots into one
        subcortexPlots = imresize(subcortexPlots,size(cortexPlot,2)/size(subcortexPlots,2)); % make same width as surface plots     
        img = [cortexPlot;subcortexPlots];
        fname_save = fullfile(folder_surfplot,[fname_surfplot,'Subcortex.png']);        

        for f = subcortexFiles
             delete(char(f));
        end

    end

    FIGURE_DISPLAY('on');