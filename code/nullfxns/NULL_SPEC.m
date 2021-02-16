function [component_design,null_spec] = NULL_SPEC(component_design_in)

% this function just splits the component_design string into a 
% component design specification and specification of null model

% INPUTS:
% component_design: string formatted as [component_design,'_',null_spec]
%
% OUTPUTS:
% component_design: 1st substring from above format
% null_spec: 2nd substring from above format

null_spec = extractAfter_(component_design_in,'_');
component_design = extractBefore_(component_design_in,'_');
% if there is no null model, then no '_' in component_design_in
% so just return [] for nullspec and input for component_design
if isempty(null_spec) & isempty(component_design)
	null_spec = [];
	component_design = component_design_in;
end