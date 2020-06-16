function txt = LABELROUND2(x,n)

% round real number x to 2 significant digits and convert to string
% used for making plot titles or overlays with stats
if ~exist('n','var')
    n=2;
end
txt = num2str(round(x,n,'significant'));
