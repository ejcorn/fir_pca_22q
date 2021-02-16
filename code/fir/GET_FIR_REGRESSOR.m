function regressor = GET_FIR_REGRESSOR(stimtiming,fin,TR,T,st)

% INPUTS:
% stimtiming: vector of stimulus onset times in seconds
% fin: number of time points after stimulus to model in the finite impulse response
% TR: scanning repetition time in seconds
% T: number of time points in scan
% st: number of time points after stimulus to start modeling
%
% OUTPUTS:
% regressor: Txfin matrix of indicators for times relative to stimulus

if ~exist('st','var')
	st = 1; % number of time points after stimulus to start modeling ... 1 in papers
end

fir = eye(fin); % finite response is a binary indicator for each time point post-stimulus
% finite response is a binary indicator for each time point post-stimulus, *start at t=1 after stimulus*

% if you have a stimulus whose predicted finite response extends beyond scan
reg_len = max(T,(max(stimtiming)/TR)+1+st+fin);
regressor = zeros(reg_len,fin); % extend regressor temporarily, then truncate it later
% otherwise make it the number of TRs in scan

for stim = 1:length(stimtiming) % iterate through stimuli
    % stimulus starting at 0s occurs during TR 1, stim at 3s occurs during
    % TR 2, etc. so use 1 + below, plus st variable for additional offset
    stimOnset = 1+st+stimtiming(stim)/TR; % get onset of response to ith stimulus in TRs, offset by st
    stimOffset = stimOnset+fin-1; % get index of where FIR ends for each stimulus, offset by st
    % add regressors together because stimuli overlap, though they will never have
    % two 1's in the same place
    regressor(stimOnset:stimOffset,:) = regressor(stimOnset:stimOffset,:) + fir;
end
regressor = regressor(1:T,:); % truncate regressor at scan duration