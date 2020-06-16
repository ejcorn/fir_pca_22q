function [scores_proj,explained_proj] = PROJECT_DATA_COEFF(X,coeff,mu)

% INPUTS:
% X: Nxp data matrix
% coeff: pxp principal components obtained by decomposing some other Mxp matrix
% with pca()
% mu: optional 1xp vector of means for the columns of X
%
% OUTPUTS:
% scores_proj: Nxp matrix X, with each column containing the projection of
% X onto the columns in coeff (each column in coeff is a PC)
% explained_proj: variance in X explained by each component of coeff
% this function will ignore any NaNs in X, and reconstruct non NaN elements
% in place

if ~exist('mu','var')
    mu = nanmean(X,1);
end

mu = repmat(mu,[size(X,1) 1]);
scores_proj = nan(size(X,1),size(coeff,2)); % initialize scores matrix
scores_proj = (X - mu) * coeff; % project centered X onto coefficients
% compute variance explained by each component in non-missing parts of X
explained_proj = [100*nanvar(scores_proj,1)/sum(nanvar(scores_proj,1))]'; 