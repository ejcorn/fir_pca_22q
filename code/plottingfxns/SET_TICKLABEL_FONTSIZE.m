function [a] = SET_TICKLABEL_FONTSIZE(xy,sz)

% INPUTS:
% xy: 'x' or 'y', corresponding to which axis tick labels you want to
% change
% sz: font size for xtick label

% OUPUTS:
% a: axis handle

% set x or y tick label font size for current subplot

if strcmp(xy,'x')
    %{
    a = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',a,'FontSize',sz);
    set(gca,'XTickLabelMode','auto');
    %}
    xl = get(gca,'XLabel');
    xlFontSize = get(xl,'FontSize');
    xAX = get(gca,'XAxis');
    set(xAX,'FontSize', sz)
    set(xl, 'FontSize', xlFontSize);
elseif strcmp(xy,'y')
    %{
    a = get(gca,'YTickLabel');  
    set(gca,'YTickLabel',a,'FontSize',sz);
    set(gca,'YTickLabelMode','auto');
    %}
    xl = get(gca,'YLabel');
    xlFontSize = get(xl,'FontSize');
    xAX = get(gca,'YAxis');
    set(xAX,'FontSize', sz)
    set(xl, 'FontSize', xlFontSize);
elseif strcmp(xy,'xy')
    xl = get(gca,'YLabel');
    xlFontSize = get(xl,'FontSize');
    xAX = get(gca,'YAxis');
    set(xAX,'FontSize', sz)
    set(xl, 'FontSize', xlFontSize);
    
    xl = get(gca,'XLabel');
    xlFontSize = get(xl,'FontSize');
    xAX = get(gca,'XAxis');
    set(xAX,'FontSize', sz)
    set(xl, 'FontSize', xlFontSize);
    
end