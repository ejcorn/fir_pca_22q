import os

os.environ['SUBJECTS_DIR'] = '/cbica/software/external/freesurfer/centos7/5.3.0/subjects/'
os.chdir('/cbica/home/cornblae/ecornblath/brain_states_22q/data/Sunetal2018_Fig1Maps/')

hemis = {'lh':'left','rh':'right'}
metrics = ['CT','SA']
# iterate through hemispheres and metrics and convert to fsaverage5
for hemi_key in hemis.keys():
	for metric in metrics:
		os.system('mri_surf2surf --hemi ' + hemi_key + ' --srcsubject fsaverage --srcsurfval ' + hemis[hemi_key] + '-'+metric+'-beta.mgh --trgsubject fsaverage5 --trgsurfval '+ hemis[hemi_key] +'-'+metric+'-beta-fsaverage5.mgh')