function idx = WHICH_MAX(x,dim)

% INPUTS:
% x: array or matrix
% dim: dimension to work across
%
% OUTPUTS:
% [~,idx] = max(x,[],dim) ... return second output of max function
% index of max values along each dimension

[~,idx] = max(x,[],dim);
idx(mean(isnan(x),dim)==1) = NaN; % if all nans in vector return nan