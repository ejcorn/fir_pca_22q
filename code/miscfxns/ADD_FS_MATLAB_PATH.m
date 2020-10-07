function ADD_FS_MATLAB_PATH()

[~,y] = system('source ~/.bash_profile; echo $FREESURFER_HOME');
addpath(genpath([y(1:end-1),'/matlab/'])); % add freesurfer matlab scripts