function [bool] = CONTAINS(str,x)

% INPUTS:
% str: string to query
% x: pattern/substr to search for
%
% OUTPUTS:
% bool: boolean indicating whether pattern (x) is in string (str)

bool = ~isempty(regexp(str,x));
