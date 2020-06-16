% construct finite impulse response regressor to measured evoked response to
% stimuli

addpaths;
cd(basedir);
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile('data',['TimeSeriesIndicators',name_root,'.mat']));
masterdir = fullfile('results',name_root);
savedir = fullfile(masterdir,'analyses','fir');
mkdir(savedir);
%% construct FIR regressor

TR=3; nTR = allScanTRs(1);
allstim = load(fullfile('data/task/idemo/stimulus/all_3col.txt'));
fin = 6; % finite response is 6 TRs long
regressor = GET_FIR_REGRESSOR(allstim(:,1),fin,TR,nTR);

save(fullfile(savedir,'RegressorFIR.mat'),'regressor','fin');
