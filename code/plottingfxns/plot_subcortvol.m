function [f_ant,f_lat] = plot_subcortvol(values,indices,sc_indices,nifti,cmap,cmin,cmax)
%%INPUTS

%values: values to base color scale off of, corresponding to indices

%indices: indices of regions corresponding to values. these indices must be
% present in the nifti matrix, but do not have to be all subcortical

%sc_indices: unique indices in nifti corresponding subcortical ROIs. This allows you to get
%rid of cortical ROIs, plot all subcortical ROIs, leaving the ones with
%nothing for 'values' as gray. You can leave out subcortical ROIs you don't
%want to appear as well by omitting them from this vector.
% nifti: 3d matrix containing brain image
%cmap: name of colormap to use
 
%cmin: optional hard minimum value for color scale
%cmax: optional hard maximum value for color scale (need both max and min)

%%OUTPUTS
%f_ant = figure handle for anterior view of subcortex
%f_med = figure handle for lateral view of subcortex

% you should check and make sure you don't have a right-left swap


if (~exist('cmap','var'))
    cmap = 'parula';
end

%%
[~,ord] = sort(indices);  % reorder indices/betas in ascending order of indices
indices = indices(ord);
values = values(ord);
values = reshape(values,[],1); indices = reshape(indices,[],1); % convert inputs to row vectors

%% convert Nx1 betas to a Nx3 matrix of RGB values

original_values = values; % store copy of original values...get in touch with your roots
if exist('cmin','var') & exist('cmax','var') % these let you specify hard limits for color scale
    if max(values) > cmax
        values(values > cmax) = cmax; % truncate max of values at cmax
    else
        values(end+1) = cmax;  % if cmax larger than max of values add it into color scaling
    end
    
    if min(values) < cmin
        values(values < cmin) = cmin; % truncate min of values at vmin
    else
        values = [cmin; values]; % if cmin smaller than min of values add it into color scaling
    end
end

% without cmin and cmax, color is just scaled based upon values
C = colormap(cmap); close; L = size(C,1);
cmap_beta_idx = round(interp1(linspace(min(values(:)),max(values(:)),L),1:L,values));
% be aware that this scaling method works worse on log/less linear data
%
if exist('cmin','var') & exist('cmax','var') 
    if max(original_values) < cmax % delete added values for scaling purposes
        cmap_beta_idx([length(values)],:) = [];
        values([length(values)],:) = [];
    end
    if min(original_values) > cmin
        cmap_beta_idx([1],:) = [];
        values([1],:) = [];
    end
else
    cmin = min(values); cmax = max(values); % set cmax and cmin to limits of values
    
end

sc = ismember(indices,unique(nifti));
indices = indices(sc); values = values(sc,:); %eliminate any non-subcortical ROIs passed in

%%

views = {[1 0 0],[0 1 0]};
n_views = length(views);
s = size(indices,1);
z = size(nifti); 

for j = 1:n_views
    if j == 1
        f_ant = figure;
        set(0,'CurrentFigure',f_ant);
    elseif j == 2
        f_lat = figure;
        set(0,'CurrentFigure',f_lat);
    end
    set(gca,'Visible','off');
    for i = 1:length(sc_indices)
        sc_grey_roi_idx = sc_indices(i);
        subcortnifti = zeros(z);
        ind = ismember(nifti,sc_grey_roi_idx);
        subcortnifti(ind) = i;
        surf = isosurface(subcortnifti,0); % if you have a LR swap, use fliplr(subcortnifti)
        p = patch('Faces',surf.faces,'Vertices',surf.vertices,'FaceColor',[0.75 0.75 0.75],'EdgeColor','none');
    end % gray background

    for i = 1:length(indices)
        sc_col_roi_idx = indices(i);
        subcortnifti = zeros(z);
        ind = ismember(nifti,sc_col_roi_idx);
        subcortnifti(ind) = i;
        surf = isosurface(subcortnifti,0);  % if you have a LR swap, use fliplr(subcortnifti)
        p = patch('Faces',surf.faces,'Vertices',surf.vertices,'EdgeColor','none','FaceVertexCData',values(i)*ones(length(surf.vertices),1),'FaceColor','flat');
    end %one by one, plot subcortical ROIs

    daspect([1 1 1])
    view(views{j})
    %colorbar('southoutside')
    colormap(cmap);
    caxis([cmin cmax])
    
    %middle head on = 1 0 0, left lat = 0 1 0, right lat= 0 -1 0
    material dull
    camlight
    lighting gouraud
end
