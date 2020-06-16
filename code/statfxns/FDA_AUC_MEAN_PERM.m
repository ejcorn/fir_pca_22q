function [A_obs,p] = FDA_AUC_MEAN_PERM(X,Y,nperms)

% INPUTS:
% X: NxT matrix whose column means define curve 1
% Y: MxT matrix whose column means define curve 2
% nperms: number of permutations to perform, determined by multiple
% comparisons correction and resulting precision needed to reject null hyp.
%
% OUTPUTS:
% A_obs: difference in the absolute value of area between two curves
% p: is the observed difference greater than the null difference?
% one-tailed test.

if ~exist('nperms','var')
    nperms = 10^4;
end

mean_fun = @(X,Y) sum(abs(nanmean(X,1) - nanmean(Y,1)));
A_obs = mean_fun(X,Y); % observed area between mean curves

%{
XY = [X;Y]; % concatenate X and Y to make null model
[NM,T] = size(XY);
samp = floor(NM/2);
A_perm = nan(nperms,1);
for perm = 1:nperms
    grp1 = randperm(NM,samp); % randomly selected half of observations from pool of both groups
    grp2 = ~ismember(1:NM,grp1); % select other half
    A_perm(perm) = mean_fun(XY(grp1,:),XY(grp2,:)); % compare areas between means of those curves
end
%}
%f=figure; plot(nanmean(X)); hold on; plot(nanmean(Y));
[N,T] = size(X);
samp = floor(N/2);
A_perm = nan(nperms,1);
for perm = 1:nperms
    X_tmp = X; Y_tmp = Y; % make new versions of X and Y
    swap = randperm(N,samp); % pick random observations to swap
    X_tmp(swap,:) = Y_tmp(swap,:); % put Y observations into new X
    Y_tmp(swap,:) = X(swap,:); % put old X observations into new Y in same positions as swap
    %f=figure; plot(nanmean(X_tmp)); hold on; plot(nanmean(Y_tmp));
    A_perm(perm) = mean_fun(X_tmp,Y_tmp); % compare areas between means of those curves
end

p = mean(A_perm>A_obs);