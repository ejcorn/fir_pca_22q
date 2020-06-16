function [data] = SCHAEFER_SURF_TO_NODE_FAST(lh_surf_data,rh_surf_data,nparc,subject)

% INPUTS:
% lh_surf_data: vertices for left hemisphere
% rh_surf_data: vertices for right hemisphere
% nparc: number of parcels
% subject: string specifying freesurfer subject, currently only fsaverage or fsaverage5

% OUTPUTS:
% data: Nxk matrix for k sets of values for N nodes in a brain parcellation, by averaging vertex values in each parcel

if ~exist('subject','var')
	subject='fsaverage5';
end

[~,k]=size(lh_surf_data);
addpath(genpath('/Applications/freesurfer/matlab'));
[annot.Rv,annot.RL, annot.Rct] = read_annotation(sprintf('data/annot/%s/rh.Schaefer2018_%dParcels_7Networks_order.annot',subject,nparc));
[annot.Lv,annot.LL, annot.Lct] = read_annotation(sprintf('data/annot/%s/lh.Schaefer2018_%dParcels_7Networks_order.annot',subject,nparc));

data = nan(nparc,k);

for j = 2:length(annot.Lct.table) % ignore medial wall
	data(j-1,:) = nanmean(lh_surf_data(find(annot.LL == annot.Lct.table(j,5)),:),1);
	data(j-1+nparc/2,:) = nanmean(rh_surf_data(find(annot.RL == annot.Rct.table(j,5)),:),1);
end