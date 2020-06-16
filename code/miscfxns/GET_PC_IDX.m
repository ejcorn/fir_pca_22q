function comp_idx = GET_PC_IDX(ncomps,XCP_folder)

% INPUTS:
% ncomps: number of components to study
% XCP_folder: name of XCP_folder specifying pipeline
%
% OUTPUTS:
% ncomps-by-1 or (ncomps+1)-by-1 cell array of chars
% e.g. {'PC1','PC2',..}
% if using a pipeline with no GSR the first component will be the global
% signal so start with PC0

if XCP_folder == 'xcp_6p_noFilter'
    comp_idx = num2cell(0:(ncomps-1));
else
    comp_idx = num2cell(1:ncomps);
end