function X = THRESHOLD(X,zdim)

% INPUTS:
% X: numeric matrix
% zdim: variable that decides absolute value threshold to apply to X

% OUTPUTS:
% X: thresholded matrix with subthreshold elements set to 0

if zdim == 3
    thresh = 0.5;
    X(abs(X) < thresh) = 0; % remove low amplitude fluctuations prior to clustering
end