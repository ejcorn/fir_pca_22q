addpaths;
cd(basedir);
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile('data',['TimeSeriesIndicators',name_root,'.mat']));
masterdir = fullfile('results',name_root);
savedir_base = fullfile(masterdir,'analyses','t1');
mkdir(savedir_base);

%% load rotations of fsaverage5 surface
subject = 'fsaverage5';
addpath(genpath('/Applications/freesurfer/matlab'))
% spin a bunch of times and save the rotation
load(fullfile(savedir_base,['SpinTestRotations',subject,'.mat']),'rotateRight','rotateLeft');

nperms = 500;

%% load reference maps -- T1 structural data from Sun et al.
% load betas
t1_metric_names = {'CT','SA'};
t1_metric_data = struct();
for met = t1_metric_names
    for hemi = {'left','right'}
        mgh_data = load_mgh(['data/Sunetal2018_Fig1Maps/',char(hemi),'-',char(met),'-beta-fsaverage5.mgh']);
        t1_metric_data.(char(met)).(char(hemi)) = squeeze(mgh_data(:,:,:,2)); % select 2nd column which is the betas for *control-22q*
    end
end
%% apply rotations to t1 maps

% apply rotations
for met = t1_metric_names % iterate through metrics and rotate surfaces
    lh_surf_data = t1_metric_data.(char(met)).left;
    rh_surf_data = t1_metric_data.(char(met)).right;
    lh_surf_data_perm = nan([size(lh_surf_data,1) nperms]);
    rh_surf_data_perm = nan([size(rh_surf_data,1) nperms]);
    testMaps_perm = nan([nparc nperms]);

    for p = 1:nperms
        fprintf('%s: Permutation %d\n',char(met),p)
        lh_surf_data_perm(:,p) = lh_surf_data(rotateLeft{p},:); % rotate all data simultaneously with each stored rotation
        rh_surf_data_perm(:,p) = rh_surf_data(rotateRight{p},:); % rotate all data simultaneously with each stored rotation
        testMaps_perm(:,p) = SCHAEFER_SURF_TO_NODE_FAST(lh_surf_data_perm(:,p),rh_surf_data_perm(:,p),nparc,subject);
    end
    t1_metric_data.(char(met)).perm = testMaps_perm;
end
%% convert spins to node data

[ct_data] = SCHAEFER_SURF_TO_NODE_FAST(t1_metric_data.CT.left,t1_metric_data.CT.right,nparc,subject);
[sa_data] = SCHAEFER_SURF_TO_NODE_FAST(t1_metric_data.SA.left,t1_metric_data.SA.right,nparc,subject);
ct_data_perm = t1_metric_data.CT.perm;
sa_data_perm = t1_metric_data.SA.perm;

save(fullfile(savedir_base,['CTSAPerm',atlasName,num2str(atlasScale),'.mat']),'ct_data_perm','sa_data_perm','ct_data','sa_data');
