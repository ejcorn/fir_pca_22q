function p = NP_TWOTAIL(dist,test_val)
% INPUTS:
% test_val: individual value being compared to distribution
% dist: Nxk matrix,distribution of values under some null model, or otherwise
%
% compute 2-tailed p-value for test value occurring in distribution
%
% OUTPUTS:
% p: kx1 vector of two-tailed p-values

p = 2*min([mean(dist <= test_val,1);mean(dist >= test_val)]);