function f = FIGURE_SIZE_CM(f,w,h)

% INPUTS:
% f: figure handle
% w: width in centimeters
% h: height in centimeters

% OUTPUTS:
% f: handle to resized figure

f.PaperUnits = 'centimeters';
f.PaperPosition = [0 0 w h];
f.PaperSize = [w h];