function [r_mat,shuffleIdx] = COLUMN_MATCH(A,B,dim,fun)

% INPUTS:
% A and B: NxK matrices
% dim: dimension. 2: match to A, 1: match to B.
% fun: anonymous function of A and B computing pairwise similarity matrix,
% e.g. corr() or @(x,y) corr(x,y).^2
%
% OUTPUTS:
% r_mat: reordered similarity matrix between columns
% shuffleIdx: indices to reorder original similarity matrix
[~,K] = size(A);

ord = 1:K;
r = fun(A,B);
[~,shuffleIdx] = max(r,[],2); % reorder
[n, bin] = histc(shuffleIdx, unique(shuffleIdx));
multiple = find(n > 1);
shuffleIdx(ismember(shuffleIdx,multiple)) = ord(~ismember(ord,shuffleIdx(~ismember(shuffleIdx,multiple))));

if dim == 1
    r_mat = fun(A(:,shuffleIdx),B);
elseif dim == 2
    r_mat = fun(A,B(:,shuffleIdx));
end