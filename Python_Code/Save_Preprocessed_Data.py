# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 11:46:29 2023

@author: kevin.tyner
"""

"""
This code does all the preprocessing of the individual MEG files, and then saves the preprocessed
file in the raw data folder.  This script allows one to analyze data once, then the remainder of
the analysis pipeline (filtering, extracting epochs, creating the source space, creating the fwd
solution, creating the inv solution, and creating the stc, followed by saving all relevant data)
is done with minimal intervention.  The preprocessed data is saved with "_proc-preproc_meg.fif" at 
the end.
"""

#%% Import modules
import os
import mne
import sys
import matplotlib.pyplot as plt
plt.ion()

#%% Set DPI to save figures
my_dpi = 96

#%% Check string for file of interest
directory = 'U:/shared/database/meg-ieeg-UNMC/'
data_read = 'rawdata'
data_write = 'derivatives'
subject_id = 'sub-unmcXXXX' 
session = 'ses-meg01'
modality = 'meg'
task = 'MEF'
side = 'ul'
run = '01'

if(side == 'None'):
    side = ''

task_id = task + side
raw_file = (directory + data_read + '/' + subject_id + '/' + session + '/' + modality + '/' + subject_id + '_' + session +
            '_task-' + task + side + '_run-' + run + '_proc-tsss_meg.fif')

if(os.path.exists(raw_file)):
    meg_path = directory + data_write + '/' + subject_id + '/' + session + '/' + task + '_MEG_output' + '_KT' ## change path depending on user
    if not(os.path.exists(meg_path)):
        os.mkdir(meg_path)
        print('MEG path made')
    else:
        #shutil.rmtree(path)
        #os.mkdir(path)
        print('MEG path exists')
    
if not(os.path.exists(raw_file)):
    print('file does not exist')
    sys.exit()
    
#%% Load file
raw = mne.io.read_raw_fif(raw_file,preload=True)

#%% Pick types
raw_meg = raw.copy().pick_types(meg=True,eog=True,ecg=True,stim=True)

#%% Set subjects directory
subjects_dir = (directory + data_write + '/' + subject_id + '/' + session + '/') 
subject = 'sample'

#%% Make transformation files
#mne.gui.coregistration(subject='sample',subjects_dir=subjects_dir,inst=raw_file)

#%% Notch filter
freqs = (60,120,180,240,300,360)
raw_meg.notch_filter(freqs=freqs,phase='zero-double')

#%% ECG SSP
if('ECG061' in raw.ch_names):
    ecg_proj,ecg_events = mne.preprocessing.compute_proj_ecg(raw,n_grad=1,n_mag=1,n_eeg=0,reject=None,ch_name='ECG061')
    
#%% EOG SSP
if('EOG062' in raw.ch_names):
    eog_proj1,eog_events1 = mne.preprocessing.compute_proj_eog(raw,n_grad=1,n_mag=1,n_eeg=0,reject=None,ch_name='EOG062')
    
#%% EOG SSP
if('EOG063' in raw.ch_names):
    eog_proj2,eog_events2 = mne.preprocessing.compute_proj_eog(raw,n_grad=1,n_mag=1,n_eeg=0,reject=None,ch_name='EOG063')
    
#%% Add SSPs
if 'ecg_proj' in locals():
    raw_meg.add_proj(ecg_proj)
if 'eog_proj1' in locals():
    raw_meg.add_proj(eog_proj1)
if 'eog_proj2' in locals():
    raw_meg.add_proj(eog_proj2)
    
#%% Apply SSPs
raw_meg.apply_proj()

#%% Copy
meg_preprocessed = raw_meg.copy()

#%% Filter
raw_meg.filter(l_freq=1,h_freq=None,phase='zero-double')

#%% ICA
ica_meg = mne.preprocessing.ICA(n_components=None,random_state=97,method='picard',max_iter=1000)
ica_meg.fit(raw_meg,picks='meg')

#%% Determine MEG components to exclude
ica_meg.plot_sources(raw_meg,start=0,stop=10,show_scrollbars=True)
ica_meg.plot_components(inst=raw_meg)

#%% Exclude MEG components
ica_meg.exclude = [] # add excluded components here
ica_rank_meg = len(ica_meg.exclude)
ica_meg.apply(inst=meg_preprocessed)

#%% Copy data
meg_examine = meg_preprocessed.copy()

#%% Filter and set epoch parameters
if(task == 'SEF'):
    meg_examine.filter(l_freq=2,h_freq=120,phase='zero-double') # change to meg_examine to look at evoked response
    tmin = -0.05
    tmax = 0.25
    baseline = (None, 0.)
elif(task == 'MEF'):
    meg_examine.filter(l_freq=0.1,h_freq=40,phase='zero-double')
    tmin = -1.0
    tmax = 0.5
    baseline = (None, -0.68)
elif(task.__contains__('lan')):
    meg_examine.filter(l_freq=1,h_freq=50,phase='zero-double')
    tmin = -1.0
    tmax = 1.0
    baseline = (-0.5, 0.)
    
#%% Create epochs
stim_channel = 'STI013'
trigger = 5
reject_dict = dict(mag=4000e-15,grad=3000e-13)

event_dict = {task_id:trigger}
events = mne.find_events(meg_examine,stim_channel=stim_channel,consecutive=False)

meg_epochs = mne.Epochs(meg_examine,events,event_id=event_dict,tmin=tmin,tmax=tmax,baseline=baseline,reject=reject_dict,
                    preload=True,reject_by_annotation=True)

#%% Create evoked MEG
meg_evoked = meg_epochs.average()
meg_evoked_topo = mne.viz.plot_evoked_topo(meg_evoked)
meg_evoked_plot = meg_evoked.plot(picks='meg',spatial_colors=True,gfp=True)

#%% Save preprocessed data
save_preprocessed_file = (directory + data_read + '/' + subject_id + '/' + session + '/' + modality + '/' + subject_id + '_' + session +
            '_task-' + task + side + '_run-' + run + '_proc-preproc_meg.fif')
meg_preprocessed.save(save_preprocessed_file,picks='all',tmin=0,tmax=None,fmt='double',overwrite=True)
