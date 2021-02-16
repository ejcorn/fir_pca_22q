function [nullTS] = UPR_BOLD_NULL(concTS,subjInd)

% INPUTS:
% concTS: (T*N)xp BOLD time series, consisting of T time points for 
% N subjects and p brain regions
% subjInd: 
%
% OUTPUTS:
% nullTS: (T*N)xp BOLD time series with each subject's time series phase randomized
% uniformly across regions

nullTS = nan(size(concTS));
for N = unique(subjInd)'
	nullTS(subjInd == N,:) = linsurr(concTS(subjInd == N,:));
end

