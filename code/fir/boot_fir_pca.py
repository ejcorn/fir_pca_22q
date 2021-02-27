import numpy as np
from sklearn.decomposition import PCA
import os
from os.path import join as opj
import scipy.io as sio
import sys
import pandas as pd

print('loaded packages')

basedir = str(sys.argv[1])
name_root = str(sys.argv[2])
component_design = str(sys.argv[3])
grp = str(sys.argv[4])
rep = int(sys.argv[5])
fin = int(sys.argv[6])
st = int(sys.argv[7])

print('got input variables')

# basedir = '/Users/Eli/Dropbox/Cornblath_Bassett_Projects/BrainStates22q/brain_states_22q'
# basedir = '/cbica/home/cornblae/ecornblath/brain_states_22q/'
# name_root = 'ScanIDSchaefer200Z1_22q_PreQC_XCP36pdespike_us_100reps'
# component_design = 'ThreatNonthreatAllStimuliStratified'
# grp = 'AllSubjects'
# rep = 1
# fin = 6
# st = 0
# os.environ['SGE_TASK_ID'] = '1'

seed = (rep-1)*100 + int(os.environ['SGE_TASK_ID']) # part of array job
print('set random seed/rep number to ' + str(seed))

#tracemalloc.start()
# load regressor and BOLD data
os.chdir(basedir)
masterdir = opj('results',name_root);
savedir_base = opj(masterdir,'analyses','fir','cpc_timecourse_fin'+str(fin)+'st'+str(st),component_design,'pncvs22qcoeff')
savedir_out = opj(savedir_base,'bootcoeff_'+grp)
os.system('mkdir -p ' + savedir_out)
m_reg = sio.loadmat(opj(savedir_base,grp+'BOLD'+component_design+'Regressor.mat'))
print('loaded data')

X = m_reg['X']
concTS = m_reg['concTS']
ncomps = int(m_reg['ncomps'])

# bootstrap sample data and recompute eigenvectors. this is detailed as a superior method in
# Peres-Neto 2003, Ecology 84, 9, 9 https://doi.org/10.1890/00-0634

np.random.seed(seed)
samp = np.random.randint(0,high=X.shape[0],size=X.shape[0]) 
X = X[samp,:]
concTS = concTS[samp,:]

# see FIR_PCA.m ... excluding regressors without responses 
# (i.e. no incorrect b/c all correct) and BOLD data for subjects 
# with no response data (i.e. actually missing)

nanmask_y = ~np.isnan(X.sum(axis=0)) 
#nanmask_x = X[:,nanmask_y.sum(axis=1)>1 # if BOLD signal is not captured by regressors at all, exclude it
nanmask_x = X[:,np.where(nanmask_y)[0]].sum(axis=1)>1 # if BOLD signal is not captured by regressors at all, exclude it
X = X[np.where(nanmask_x)[0],:] # delete rows that can't be modeled in this sample

# find duplicate columns introduced and remove
df = pd.DataFrame(X.T) 
dup_mask = df.duplicated(keep=False).values
X[:,np.where(dup_mask)[0]] = np.nan
nanmask_y = ~np.isnan(X.sum(axis=0)) 

X = X[:,np.where(nanmask_y)[0]] # delete columns that can't be modeled in this sample
concTS = concTS[np.where(nanmask_x)[0],:] # delete rows (BOLD) that can't be modeled

print('prepped BOLD and regressor')
# regress concTS onto X using np.linalg.lstsq(a, b, rcond=-1): 
# a x = b by computing a vector x that minimizes the Euclidean 2-norm || b - a x ||^2.

betaFIR,residuals,rank,s = np.linalg.lstsq(X,concTS,rcond=None)
#TS_reconstruct =  np.matmul(X,betaFIR) # reconstruct TS as predicted values of this regression

print('performed regression')

# decompose task-related variance in TS_reconstruct
pca = PCA(n_components=ncomps)
pca.fit(np.matmul(X,betaFIR)) # reconstruct TS as predicted values of this regression then decompose

print('performed PCA')

sio.savemat(file_name=opj(savedir_out,'BootstrapPCA_Rep'+str(seed)+'.mat'),mdict={'coeff':pca.components_.T,'samp':samp,'nanmask_y':nanmask_y,'nanmask_x':nanmask_x,'explained':pca.explained_variance_,'rank':rank})
#snapshot = tracemalloc.take_snapshot()

print('saved data')