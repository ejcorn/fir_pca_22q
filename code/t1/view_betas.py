import nibabel as nib
import numpy as np
import os
from os.path import join as opj
import scipy.io as sio
from glob import glob
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib import rc
from matplotlib.colors import ListedColormap
from matplotlib.colors import Normalize
import mayavi as my
import surfer
from scipy import ma
from scipy import interp

import sys
sys.path.append('code/visualize/ejcbrain')
import ejcbrain as EB

hemi = 'lh'
views = ['lat','med']
thrsh = -10
subject_id = 'fsaverage5'
surf = 'pial'
cbar = False
cmap_array = np.load('data/colors/ejc_custom2.npy')
cmap = ListedColormap(cmap_array)

font = {'family' : 'arial',
        'size'   : 6}
rc('font', **font)

# create a class to make an asymmetric colorbar
class MidpointNormalize(Normalize):
    def __init__(self, vmin, vmax, midpoint=0, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        normalized_min = max(0, 1 / 2 * (1 - abs((self.midpoint - self.vmin) / (self.midpoint - self.vmax))))
        normalized_max = min(1, 1 / 2 * (1 + abs((self.vmax - self.midpoint) / (self.midpoint - self.vmin))))
        normalized_mid = 0.5
        x, y = [self.vmin, self.midpoint, self.vmax], [normalized_min, normalized_mid, normalized_max]
        return ma.masked_array(interp(value, x, y))

for m in ['CT','SA']:
	for view in views:
		b = nib.freesurfer.mghformat.load('data/Sunetal2018_Fig1Maps/left-'+m+'-beta.mgh') 
		p = nib.freesurfer.mghformat.load('data/Sunetal2018_Fig1Maps/left-'+m+'-sig.mgh') 
		# p = nib.freesurfer.mghformat.load('data/Sunetal2018_Fig1Maps/left-CT-sig.mgh') 
		b = np.squeeze(np.asarray(b.get_data()[:,:,:,1])) # select 22q-control predictor betas
		p = np.squeeze(np.asarray(p.get_data()))

		data_plot = -b
		clim = np.abs(b).max()
		cmax = b.max()
		cmin = b.min()
		fig = my.mlab.figure(size=(340,340))
		fig = my.mlab.gcf()
		brain = surfer.Brain(subject_id, hemi, surf,figure=fig,views=view,background='white', alpha = 1)
		brain.add_data(data_plot,min=cmin,max=cmax,mid=0,thresh = thrsh, colormap= cmap,alpha=1,colorbar=cbar)
		my.mlab.savefig(figure=fig,filename=opj('data','Sunetal2018_Fig1Maps',view+m+subject_id+'Betas.png'),magnification=1)

	x = np.vstack([np.linspace(cmin,cmax,100) for y in np.arange(3)]) 
	fig, ax = plt.subplots()
	cax=ax.imshow(x,cmap=cmap,norm=MidpointNormalize(vmin=cmin,vmax=cmax,midpoint=0))
	fig.colorbar(cax,ticks=[cmin,0,cmax],orientation='horizontal')
	fig.set_size_inches(1.2,2)
	fig.savefig(opj('data','Sunetal2018_Fig1Maps',view+m+subject_id+'BetaColormap.pdf'))
