function mav = MAX_ABS_VALUE(x,dim)
% INPUTS:
% x: matrix
% dim: dimension
% 
% OUTPUT:
% mav: along specified dimension, return the value that is the maximum
% absolute value (i.e. in each row or column)

if dim == 2
    ind = sub2ind(size(x),[1:size(x,1)]',WHICH_MAX(abs(x),2));
elseif dim == 1
    ind = sub2ind(size(x),WHICH_MAX(abs(x),1)',[1:size(x,2)]');
end
mav = nan(size(ind,1),1);
mav(~isnan(ind)) = x(ind(~isnan(ind))); % replace non-NaN values in place
    