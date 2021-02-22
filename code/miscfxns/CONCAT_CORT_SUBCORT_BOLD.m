function [concTS_all,nparc_all,cort_indices,sc_indices_combined] = CONCAT_CORT_SUBCORT_BOLD(concTS,atlasNameSubcortex,atlasScaleSubcortex,atlasScale,name_root)

% INPUTS:
% concTS: concatenated BOLD time series matrix that is (subjects*# TRs per scan)x(regions or parcels)
% atlasNameSubcortex: character containing name of subcortical atlas. here, we use 'HarvardOxford'
% atlasScaleSubcortex: double/int containing parcellation scale of subcortical atlas
% atlasScale: cortical atlas scale, number of parcels
% name_root: name_root variable indexing project master folder
%
% OUTPUTS:
% concTS_all: concatenated BOLD time series with subcortical regions added
% nparc_all: new number of regions in concTS
% cort_indices: indices of cortical parcels, to be applied to columns of concTS_all
% sc_indices_combined: indices of subcortical parcels, to be applied to columns of concTS_all
	
nparc = size(concTS,2);
[nifti,sc_indices] = RETURN_NII_SUBCORT(atlasNameSubcortex,atlasScaleSubcortex);
extralab = extractAfter_(name_root,num2str(atlasScale));
% load time series data from
concTS_subcort = load(fullfile('data',sprintf(['ConcTimeSeriesCPCA_ID%s%d',extralab,'.mat'],atlasNameSubcortex,atlasScaleSubcortex)));
%% concatenate harvard oxford subcortex with schaefer subcortex
concTS_all = [concTS concTS_subcort.concTS(:,sc_indices)];
cort_indices = 1:nparc; % indices for schaefer cortical parcels
nparc_all = size(concTS_all,2);
sc_indices_combined = (nparc+1):nparc_all; % location of subcortical nodes in concatenated matrix