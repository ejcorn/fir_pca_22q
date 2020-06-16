addpaths;
cd(basedir);
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile('data',['TimeSeriesIndicators',name_root,'.mat']));
masterdir = fullfile('results',name_root);
savedir = fullfile(masterdir,'analyses','t1');
mkdir(savedir);

%% load separate cluster solutions

cluster_dir = fullfile(masterdir,'analyses','centroids','pncvs22qseparatecluster');
PNC.(['k',num2str(numClusters)]) = load(fullfile(cluster_dir,['PNCOnlyClusterCentroids_k',num2str(numClusters),name_root,'.mat']));
q22.(['k',num2str(numClusters)]) = load(fullfile(cluster_dir,['22qOnlyClusterCentroids_k',num2str(numClusters),name_root,'.mat']));

% compute difference in matched centroid maps between the two clustering solutions
q22_centroids = q22.(['k',num2str(numClusters)]).kClusterCentroids;
PNC_centroids = PNC.(['k',num2str(numClusters)]).kClusterCentroids;
q22MinusPNC_CentroidDifference = q22_centroids - PNC_centroids;

%% compute rotations of fsaverage5 surface
subject = 'fsaverage5';

% read in fsaverage5 surface
[Lvertices, Lfaces] = freesurfer_read_surf(sprintf('data/surf/%s/lh.pial',subject));
[Rvertices, Rfaces] = freesurfer_read_surf(sprintf('data/surf/%s/rh.pial',subject));

% spin a bunch of times and save the rotation
parpool(4);

nperms = 500;
rotateRight = cell(nperms,1);
rotateLeft = cell(nperms,1);
parfor p = 1:nperms
    sprintf('Permutation %d',p) % repeat 1 permutation nperms times and store rotations (doing this to paralellize
    [~,~,Ir_i,Il_i] = SpinPermuFS_EJCFAST(Lvertices,Rvertices,1);
    rotateRight{p} = Ir_i{1};
    rotateLeft{p} = Il_i{1};
end
save(fullfile(savedir,['SpinTestRotations',subject,'.mat']),'rotateRight','rotateLeft');

%% apply rotations to data
testMaps = [q22_centroids PNC_centroids];
% convert centroid difference maps to fsaverage5 surface
[lh_surf_data,rh_surf_data] = SCHAEFER_NODE_TO_SURF_FAST(testMaps,subject);

% apply rotations
lh_surf_data_perm = nan([size(lh_surf_data) nperms]);
rh_surf_data_perm = nan([size(rh_surf_data) nperms]);
node_data_perm = nan([size(testMaps) nperms]);

for p = 1:nperms
    sprintf('Permutation %d',p)
    lh_surf_data_perm(:,:,p) = lh_surf_data(rotateLeft{p},:); % rotate all data simultaneously with each stored rotation
    rh_surf_data_perm(:,:,p) = rh_surf_data(rotateRight{p},:); % rotate all data simultaneously with each stored rotation
    node_data_perm(:,:,p) = SCHAEFER_SURF_TO_NODE_FAST(lh_surf_data_perm(:,:,p),rh_surf_data_perm(:,:,p),nparc,subject);
end

%% load reference maps -- T1 structural data from Sun et al.
d = dir('data/Sunetal2018_Fig1Maps/*beta*.mgh');
fnames = strcat({d.folder}','/',{d.name}');
mgh_data = cellfun(@(x) load_mgh(x),fnames,'Uni',false);
mgh_data = cellfun(@(x) squeeze(x(:,:,:,2)),mgh_data,'UniformOutput',false);
t1_data = [mgh_data{:}];
[ct_data] = SCHAEFER_SURF_TO_NODE_FAST(t1_data(:,1),t1_data(:,3),nparc,'fsaverage');
[sa_data] = SCHAEFER_SURF_TO_NODE_FAST(t1_data(:,2),t1_data(:,4),nparc,'fsaverage');

%% null vs actual variance explained

refmap_cell = {ct_data,sa_data};
r2_actual = cellfun(@(x) corr(x,testMaps).^2,refmap_cell,'UniformOutput',false);
r2_null = cell(length(refmap_cell),1);
p_spin = cell(length(refmap_cell),1);
for j = 1:length(refmap_cell)
    r2_null{j} = nan(nperms,size(testMaps,2));
    for p = 1:nperms
        r2_null{j}(p,:) = corr(refmap_cell{j},node_data_perm(:,:,p),'rows','pairwise').^2;
    end
    % p value interpretation: probability that you observe greater R^2 in
    % null
    p_spin{j} = mean(r2_actual{j} < r2_null{j},1); 
end
p_spin{2}
%% plot
sqrt(vertcat(r2_actual{:}))
vertcat(p_spin{:})
f=figure;
plot(refmap_cell{2},testMaps(:,5),'b.');
