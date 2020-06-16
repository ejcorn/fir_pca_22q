import nibabel as nib
import numpy as np
import os
from os.path import join as opj
import scipy.io as sio
from glob import glob

class node_to_vertex:

	def __init__(self,annot_path,node_data):
		self.node_data = node_data
		if len(node_data.shape) < 2:
			self.node_data = np.expand_dims(node_data,axis=1) # make node-by-condition, even if 1 condition
		self.annot_path = annot_path

	def getvertdata_lausanne(self,hemi = ''):

		node_data = self.node_data
		nparc = node_data.shape[0] # get number of ROIs

		# assign lausanne "scale" parameter based on number of ROIs, which determines name of annot file
		# parc_ind specifies where in the MATLAB cell holding region names that a specific scale can be found
		if nparc < 100:
			scl = '36'
			parc_ind = 0
		elif nparc > 100 and nparc < 200:
			scl = '60'
			parc_ind = 1
		elif nparc > 200 and nparc < 300:
			scl = '125'
			parc_ind = 2
		elif nparc > 400 and nparc < 500:
			scl = '250'
			parc_ind = 3

		annot_fname = opj(self.annot_path,hemi+'.myaparc_'+scl+'.annot')
		labels, ctab, names = nib.freesurfer.read_annot(annot_fname)

		# load lausanne region names, which exist in a MATLAB cell
		roinames = sio.loadmat(opj(self.annot_path,'human_regionNames.mat'))['roinames'][0][parc_ind]
		nparc = len(roinames)
		for i in np.arange(0,nparc):
			roinames[i] = str(roinames[i][0][0])

		for i in np.arange(0,len(names)):   # add hemisphere to names
			names[i] = hemi + "_" + names[i]
			
		present_label_mask = np.in1d(roinames,names)
		roinames = roinames[present_label_mask]
		missinglabels = [False]*len(names)
		idx_init = [0]*len(names)
		for i in np.arange(0,len(names)):
			idx_init[i] = np.where(roinames == names[i])[0]
			missinglabels[i] = np.size(idx_init[i]) == 0

		#labels[np.in1d(labels,missinglabels[0])] = -1000

		# get rid of [] list elements -- this is a terrible way to do this but i don't know python
		presentlabels = np.where(np.invert(missinglabels))
		missinglabels = np.where(missinglabels)
		idx = [0]*np.size(presentlabels)
		for i in np.arange(0,np.size(presentlabels)):
			idx[i] = idx_init[presentlabels[0][i]]

		idx = np.array(idx)
		vtx_data = np.double(node_data[present_label_mask,:])
		vtx_data = np.squeeze(vtx_data[idx,:])
		for i in np.arange(0,len(missinglabels[0])):
			vtx_data = np.insert(vtx_data, obj = missinglabels[0][i],values = -10000, axis = 0)

		if len(vtx_data.shape) < 2:
			vtx_data = np.expand_dims(vtx_data,axis=1) # make node-by-condition, even if 1 condition
		vtx_data = vtx_data[labels,:]
		return vtx_data

	def getvertdata_Yeo(self,hemi=''):

		# INPUTS:
		# node_data: data associated with Yeo parcels. input data for all parcels to automatically detect parcellation scale
		# assumes first half of node_data is left hemi, second half is right hemi
		# hemi: 'lh' or 'rh', which hemisphere to plot
		#
		# OUTPUTS
		# vtx_data: node_data values assigned to respective surface vertices

		node_data = self.node_data
		nparc = node_data.shape[0]
		annot_fname = opj(self.annot_path,hemi+'.Schaefer2018_'+str(nparc)+'Parcels_7Networks_order.annot')
		labels, ctab, names = nib.freesurfer.read_annot(annot_fname)
		
		if hemi == 'lh':
			node_data = node_data[np.arange(int(nparc/2)),:]
		if hemi == 'rh':
			node_data = node_data[int(nparc/2):nparc,:]
		
		# add -1000 as first value in every hemisphere's node data ( this corresponds to medial wall)
		vtx_data = np.concatenate((np.tile(np.array([-999]),(1,node_data.shape[1])),node_data))
		vtx_data = vtx_data[labels,:]
		return vtx_data

	def getvertdata_BN(self,hemi=''):

		# INPUTS:
		# node_data: data associated with brainnetome parcels. input data for all parcels to automatically detect parcellation scale
		# assumes node data alternates between left and right hemisphere for matching parcels
		# hemi: 'lh' or 'rh', which hemisphere to plot
		#
		# OUTPUTS
		# vtx_data: node_data values assigned to respective surface vertices

		node_data = self.node_data
		nparc = node_data.shape[0]
		annot_fname = opj(self.annot_path,hemi+'.BN_Atlas.annot')
		labels, ctab, names = nib.freesurfer.read_annot(annot_fname)
		
		# if hemi == 'lh':
		# 	hemi_nodes = np.where(['_L' in name for name in names])[0]
		# if hemi == 'rh':
		# 	hemi_nodes = np.where(['_R' in name for name in names])[0]
		
		# node_data = node_data[hemi_nodes,:]
		
		# add -1000 as first value in every hemisphere's node data ( this corresponds to medial wall)
		labels[labels == -1] = 0 # 0 and -1 in labels are not mapped to any cortical ROIs
		vtx_data = np.concatenate((np.tile(np.array([-1000]),(1,node_data.shape[1])),node_data))
		vtx_data = vtx_data[labels,:]
		return vtx_data

	def getvertdata(self,atlas='',hemi=''):

		# INPUTS:
		# atlas: string specifying which atlas you are using and thus which vertex mapping function to usse
		# hemi: which hemisphere to plot

		if 'Yeo' in atlas or 'Schaefer' in atlas or 'yeo' in atlas:
			return self.getvertdata_Yeo(hemi)
		elif 'Laus' in atlas:
			return self.getvertdata_lausanne(hemi)
		elif 'Brainnetome' in atlas:
			return self.getvertdata_BN(hemi)

def getannot_Yeo(nparc,hemi='',annot_path='data/annot'):
	annot_fname = opj(annot_path,hemi+'.Schaefer2018_'+str(nparc)+'Parcels_7Networks_order.annot')
	labels, ctab, names = nib.freesurfer.read_annot(annot_fname)
	return labels,ctab,names


import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import mayavi as my
import surfer
#my.mlab.options.offscreen = True

def surfplot_2xK(nodeData,hemi,atlas,ttls,savedir,clim_type='each',cbar_yn=True,cmap = 'plasma'):

	# INPUTS:
	# nodeData: NxK data for N nodes and K activation maps
	# hemi: hemisphere to plot
	# atlas: string specifiying atlas of node data
	# ttls: length K list of plot titles
	# savedir: where to output temporary files which will be deleted
	# clim_type: optional limit for color axis, symmetric about 0. 
	# 'each' (default) means each map is scaled independently
	# 'all' means set color scales based on all values in nodeData
	# a numeric input will be set as the color limit
	# cbar_yn: if True, put a colorbar on the bottom row of plots
	
	views=['lat','med']
	thrsh = -10
	subject_id = 'fsaverage'
	surf = 'pial'
	tmp = 'TEMP_BRAIN'
	
	N,K = nodeData.shape
	ntv = node_to_vertex(annot_path='data/annot/',node_data=nodeData)
	vtx_data = ntv.getvertdata(atlas,hemi=hemi)	
	img=dict()
	for act_map in range(K):
		data_plot = vtx_data[:,act_map]
		if clim_type == 'each':
			clim = np.max(np.abs(data_plot[data_plot>-999])) 
		elif clim_type == 'all':
			clim = np.max(np.abs(vtx_data[vtx_data > -999]))
		else:
			clim = clim_type
		print('color limits:' +str(clim))
		for view in views:
			cbar = False
			if view == 'med':
				cbar = cbar_yn 
			fig = my.mlab.figure(size=(340,340))
			fig = my.mlab.gcf()
			brain = surfer.Brain(subject_id, hemi, surf,figure=fig,views=view,background='white', alpha = 1)
			brain.add_data(data_plot, min = -clim, max = clim, thresh = thrsh, colormap= cmap,alpha=1,colorbar=cbar)       
			fname = tmp + str(act_map) + view + hemi + '.png'
			#img[str(act_map)+view] = my.mlab.screenshot(figure=fig, mode='rgba', antialiased=True)
			my.mlab.savefig(figure=fig,filename=opj(savedir,fname),magnification=1)
	my.mlab.close(all=True)

	fig = plt.figure(figsize = [K,2])
	arial = {'fontname':'Arial'}
	for view in np.arange(1,len(views)+1):
		for act_map in np.arange(K):
			fname = tmp + str(act_map) + views[view-1] + hemi + '.png'
			img = mpimg.imread(opj(savedir,fname))
			plt.subplot(2,K,1+act_map + K*(view-1))
			#imgplot = plt.imshow(img[str(act_map)+views[view-1]],aspect='auto')
			imgplot = plt.imshow(img,aspect='auto')
			if view == 1:
				ttl = ttls[act_map]
				plt.title(ttl,fontsize=8,fontweight='bold',**arial)
			plt.axis('off')
	plt.subplots_adjust(hspace=0, wspace=0.05)
	[os.remove(f) for f in glob(opj(savedir,tmp + '*' + hemi + '.png'))]
	return(plt.gcf())

def surfplot_2x2(nodeData,atlas,savedir,ttl='',clim_type='all',cmap='plasma'):
	# INPUTS:
	# nodeData: Nx1 data for N nodes
	# atlas: string specifiying atlas of node data
	# ttl: title to go above entire plot
	# savedir: where to output temporary files which will be deleted
	# clim_type: optional limit for color axis, symmetric about 0. 
	# 'all' (default) means set color scales based on all values in nodeData
	# a numeric input will be set as the color limit

	if clim_type == 'all':
		clim = np.round(np.max(np.abs(nodeData)),decimals=2)
		print('color limits:' +str(clim))
		if clim ==0:
			clim = 0.2
			print('setting clim to 0.2 by default')
	else:
		clim = clim_type
	hemis = ['lh','rh']
	views = ['lat','med']
	thrsh = -10
	subject_id = 'fsaverage'
	surf = 'pial'
	cbar = False

	for hemi in hemis:
		ntv = node_to_vertex(annot_path='data/annot/',node_data=nodeData)
		vtx_data = ntv.getvertdata(atlas,hemi=hemi)	
		for view in views:  
			fig = my.mlab.figure(size=(340,340))
			fig = my.mlab.gcf()
			brain = surfer.Brain(subject_id, hemi, surf,figure=fig,views=view,background='white', alpha = 1)
			brain.add_data(vtx_data, min = -clim, max = clim, thresh = thrsh, colormap=cmap, alpha=.8,colorbar=cbar,time_label=None)    
			fname = 'nodeData' + view + hemi + '.png'
			my.mlab.savefig(figure=fig,filename=opj(savedir,fname),magnification=1)	
	my.mlab.close(all=True)

	#Arrange brains into grid

	arial = {'fontname':'Arial'}

	plt.figure(figsize = [1.7,1.7])
	for H,hemi in enumerate(hemis):
		for V, view in enumerate(views):    #Arrange brains into grid
			fname = opj(savedir,'nodeData' + view + hemi + '.png')
			img = mpimg.imread(fname)
			plt.subplot(2,2,(H + 2*V+1))
			imgplot = plt.imshow(img,aspect='auto')
			plt.axis('off')

	plt.subplots_adjust(hspace=0, wspace=0.05)
	plt.suptitle(ttl,family='arial',size=8,weight='bold')
	#delete intermediate files
	[os.remove(f) for f in glob(opj(savedir,'nodeData*.png'))]
	return(plt.gcf())

def surfplot_2x2xk(nodeData,atlas,ttls,savedir,clim_type='each',cmap = 'plasma'):

	# INPUTS:
	# nodeData: NxK data for N nodes and K activation maps
	# atlas: string specifiying atlas of node data
	# ttls: length K list of plot titles
	# savedir: where to output temporary files which will be deleted
	# clim_type: optional limit for color axis, symmetric about 0. 
	# 'each' (default) means each map is scaled independently
	# 'all' means set color scales based on all values in nodeData
	# a numeric input will be set as the color limit

	# will plot 2x2 plots of both hemispheres, medial and lateral

	num_brains = nodeData.shape[1]
	for j in np.arange(num_brains): # make full brain 2x2 plots
		f=surfplot_2x2(nodeData[:,j],atlas,savedir,ttl=ttls[j],clim_type=clim_type,cmap=cmap)
		fname = 'TEMP_2x2_' + str(j)+".png"
		f.savefig(opj(savedir,fname),dpi=500,bbox_inches='tight',pad_inches=0)

	plt.figure(figsize = [1.72*num_brains,1.7])
	for j in np.arange(num_brains): # tile 2x2 plots
		fname = 'TEMP_2x2_' + str(j)+".png"
		img = mpimg.imread(opj(savedir,fname))
		plt.subplot(1,num_brains,j+1)
		imgplot = plt.imshow(img,aspect='equal')
		plt.axis('off')
	plt.subplots_adjust(hspace=0, wspace=0)

	[os.remove(opj(savedir,'TEMP_2x2_' + str(j)+".png")) for j in np.arange(num_brains)]

	return(plt.gcf())
