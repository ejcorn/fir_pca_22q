import nibabel as nib
import numpy as np
import os
from os.path import join as opj
import scipy.io as sio
from glob import glob
import neurosynth as ns

homedir = '/Users/Eli/Dropbox/Cornblath_Bassett_Projects/BrainStates22q/fir_pca_22q/'
resultsdir = 'results/CPCA_IDSchaefer200Z1xcp_6p_noFilter/analyses/fir/subject_fir_correct_incorrect_pca/cpc_timecourse/ThreatNonthreatAllStimuliStratified/pncvs22qcoeff/'
savedir = opj(homedir,resultsdir,'nifti')
os.system('mkdir -p '+savedir)
m = sio.loadmat(opj(homedir,resultsdir,'FIRGroupThreatNonthreatAllStimuliStratified_CPCAComponents.mat'))

coeffs = m['nodeDataAll']

ncomps = coeffs.shape[1]-1

nparc = 200
schaefer = nib.load('data/nifti/Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.nii.gz')
ho = nib.load('data/nifti/HarvardOxford/HarvardOxfordMNI.nii.gz')
s_img = schaefer.get_fdata() # schaefer nifti
h_img = ho.get_fdata() # harvard oxford nifti
PC = 0
for PC in np.arange(ncomps):
	n_img = np.zeros(s_img.shape) # new blank nifti
	for j in np.arange(nparc):
		n_img[s_img==(j+1)] = coeffs[j,PC] # replace each node with schaefer data from each pc
	ho_inds = np.unique(h_img)[-14:]
	for j in np.arange(ho_inds.shape[0]):
		n_img[h_img==ho_inds[j]] = coeffs[j+nparc,PC]
		

	n_img = nib.Nifti1Image(n_img, ho.affine)
	n_img.to_filename(opj(savedir,'PC'+str(PC)+'Unthresholded.nii'))

from neurosynth.base.dataset import Dataset
from neurosynth.analysis import decode
dataset = Dataset('data/database.txt')
dataset.add_features('data/features.txt')
dataset.save('data/dataset.pkl')

dataset = Dataset.load('data/dataset.pkl')
decoder = decode.Decoder(dataset)
result = decoder.decode([opj(savedir,'PC'+str(PC)+'.nii') for PC in np.arange(ncomps)],save = opj(savedir,'decoding_results.txt'))

for PC in np.arange(ncomps):
	print('PC'+str(PC))
	print(result.iloc[:,PC].sort_values().iloc[:10])
	print(result.iloc[:,PC].sort_values().iloc[-10:])