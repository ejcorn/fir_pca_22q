cd(basedir); addpath(genpath('code'));
addpaths;
cd(basedir);
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile('data',['TimeSeriesIndicators',name_root,'.mat']));
masterdir = fullfile('results',name_root);
savedir = fullfile(masterdir,'analyses','t1');
mkdir(savedir);

%% compute rotations of fsaverage5 surface
subject = 'fsaverage5';
ADD_FS_MATLAB_PATH()

% read in fsaverage5 surface
[Lvertices, Lfaces] = freesurfer_read_surf(sprintf('data/surf/%s/lh.pial',subject));
[Rvertices, Rfaces] = freesurfer_read_surf(sprintf('data/surf/%s/rh.pial',subject));

% spin a bunch of times and save the rotation
%parpool(4);

nperms = 500;
rotateRight = cell(nperms,1);
rotateLeft = cell(nperms,1);
for p = 1:nperms
    sprintf('Permutation %d',p) % repeat 1 permutation nperms times and store rotations (doing this to paralellize
    [~,~,Ir_i,Il_i] = SpinPermuFS_EJCFAST(Lvertices,Rvertices,1);
    rotateRight{p} = Ir_i{1};
    rotateLeft{p} = Il_i{1};
end
save(fullfile(savedir,['SpinTestRotations',subject,'.mat']),'rotateRight','rotateLeft');