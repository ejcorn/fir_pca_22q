
% one of these two must contain load_nii functions
niftitoolbox_path1 = '/cbica/home/cornblae/ecornblath/matlab/brainmapping2/NIfTI_20140122';
niftitoolbox_path2 = '~/Dropbox/ToolBoxes/NIfTI_20140122';

BCTPath = '/cbica/home/cornblae/ecornblath/matlab/BCT/';
NCTPath = '/cbica/home/cornblae/ecornblath/matlab/control_fc/NCT/';

addpath(genpath(BCTPath));
addpath(genpath(NCTPath));
addpath(genpath(niftitoolbox_path1)); addpath(genpath(niftitoolbox_path2));

a = clock; rng(a(6));