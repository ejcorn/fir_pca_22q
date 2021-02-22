function [f,X_subj] = PLOT_DESIGN_MATRIX(design_str,subj,subjInd)

% INPUTS:
% design_str: structure containing FIR design matrix, subject and stimulus labels
% subj: optional, select just one subject
% subjInd: if you give subj input, also give subjInd
%
% OUTPUTS:
% f: figure handle to plot showing design matrix
% X_subj: design matrix plotted

m = 100; % increase height of time line
b = 1000; % shift subject line 
s = 40; % density of tick labels

if (sum(design_str.X(:,1)) ~= size(design_str.X,1)) & (size(design_str.X,2) < length(design_str.columnLabelsSubject))
    % add intercept if it's missing (accommodates use of this function
    % around line 140 of code/fir/fir_design.m)
    design_str.X = [ones(size(design_str.X,1),1) design_str.X];
end

if exist('subj','var')
	subjInd_columns = [1 find(design_str.columnLabelsSubject==subj)]; % keep intercept too
	design_str.X = design_str.X(subjInd==subj,subjInd_columns);
	design_str.columnLabelsSubject = design_str.columnLabelsSubject(subjInd_columns);
	design_str.columnLabelsStim = design_str.columnLabelsStim(subjInd_columns);
	design_str.columnLabelsResp = design_str.columnLabelsResp(subjInd_columns);
	design_str.columnLabelsTime = design_str.columnLabelsTime(subjInd_columns);
	m=30; b = 10; s = 1;
end

[r,c] = size(design_str.X);
disp([r,c])
xt = 1:c; sel = 1:s:c;
xt = xt(sel);
f=figure;
imagesc(design_str.X);
xticks(xt); xticklabels(design_str.columnLabelsResp(sel)); xtickangle(90)
yticklabels({});
ax1=gca;
ax2 = axes('Position', get(ax1, 'Position'),'Color', 'none');
set(ax2, 'XAxisLocation', 'top');
set(ax2, 'XLim', get(ax1, 'XLim'),'YLim',get(ax1,'YLim'));
set(ax2,'XTick',xt,'XTickLabel',design_str.columnLabelsStim(sel),'XTickLabelRotation',90)
hold on; plot(1:c,m*design_str.columnLabelsTime,'r');
hold on; plot(1:c,b+(m/10)*design_str.columnLabelsSubject,'g')

X_subj = design_str.X;