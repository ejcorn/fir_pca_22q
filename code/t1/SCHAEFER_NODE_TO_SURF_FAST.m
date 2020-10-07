function [lh_surf_data,rh_surf_data] = SCHAEFER_NODE_TO_SURF_FAST(data,subject)

% INPUTS:
% data: Nxk matrix for k sets of values for N nodes in a brain parcellation

% OUTPUTS:
% lh_surf_data: vertices for left hemisphere
% rh_surf_data: vertices for right hemisphere

if ~exist('subject','var')
	subject='fsaverage5';
end

[N,k] = size(data);
ADD_FS_MATLAB_PATH()

[annot.Rv,annot.RL, annot.Rct] = read_annotation(sprintf('data/annot/%s/rh.Schaefer2018_%dParcels_7Networks_order.annot',subject,N));
[annot.Lv,annot.LL, annot.Lct] = read_annotation(sprintf('data/annot/%s/lh.Schaefer2018_%dParcels_7Networks_order.annot',subject,N));

lh_data = vertcat(repmat(NaN,[1 k]),data(1:(N/2),:));
rh_data = vertcat(repmat(NaN,[1 k]),data((N/2)+1:end,:));

lh_surf_data = nan(size(annot.Lv,1),k);
rh_surf_data = nan(size(annot.Rv,1),k);

l_labels = arrayfun(@(x) find(annot.Lct.table(:,5) == x),annot.LL);
r_labels = arrayfun(@(x) find(annot.Rct.table(:,5) == x),annot.RL);

for k_i = 1:k
    lh_surf_data(:,k_i) = lh_data(l_labels,k_i);
    rh_surf_data(:,k_i) = rh_data(r_labels,k_i);
end

% save('data/tmp.mat','lh_surf_data','rh_surf_data') % confirmed that
% surface plots look as expected using this script
