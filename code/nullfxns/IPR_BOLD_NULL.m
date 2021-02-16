function [nullTS] = IPR_BOLD_NULL(concTS,subjInd)

% INPUTS:
% concTS: (T*N)xp BOLD time series, consisting of T time points for 
% N subjects and p brain regions
% subjInd: 
%
% OUTPUTS:
% nullTS: (T*N)xp BOLD time series with each region's time series phase randomized within
% each subject

nullTS = nan(size(concTS));
for N = unique(subjInd)'
	nullTS(subjInd == N,:) = linsurr_ind(concTS(subjInd == N,:));
end

