# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 10:51:20 2023

@author: kevin.tyner
"""

#%% Import external functions
import sys
sys.path.insert(0,'path/to/functions/') # change path

#%% Import modules
import os
import os.path as op
import numpy as np
import mne
import re
import matplotlib.pyplot as plt
from scipy.io import savemat
from Preprocessed_Dict import Preprocessed_Dict
import gc
plt.ion()

#%% Set dict
my_dict = Preprocessed_Dict()

#%% Set DPI to save figures
my_dpi = 120 # originally 96

#%% Determine Windows or Linux
if(os.path.exists('U:/shared/database/meg-ieeg-UNMC/')):
    directory = 'U:/shared/database/meg-ieeg-UNMC/'
elif(os.path.exists('/mnt/udata/shared/database/meg-ieeg-UNMC/')):
    directory = '/mnt/udata/shared/database/meg-ieeg-UNMC/'
    
#%% Read from server and sort numerically
data_location = directory + 'rawdata/'
Subjects = os.listdir(data_location)

#%% Remove non-unmc patients
remove_indices = []
for idx in range(len(Subjects)):
    if not('unmc' in Subjects[idx]):
        remove_indices.append(idx)
somelist = [i for j,i in enumerate(Subjects) if j not in remove_indices] # delete patients that are not in the database
Subjects = somelist
del somelist,remove_indices

#%% Sort Subjects
Subjects.sort(key=lambda f: int(re.sub('\D','',f)))

#%% Delete subjects I cannot access
del Subjects[89:99]
del Subjects[246:] # del to end

#%% New location for subjects_dir
location = directory + 'derivatives/'

#%% Loop through patients to get all .fif files
for subject_idx in range(len(Subjects)): 
    subjects_data = data_location + Subjects[subject_idx] + '/ses-meg01'
    subjects_dir = location + Subjects[subject_idx] + '/ses-meg01'
    Patient_files = subjects_data + '/meg'
    Subject = Subjects[subject_idx]
    subject = 'sample'
    fif_files = list()
    for file in os.listdir(Patient_files):
        if file.endswith('.fif'):
            fif_files.append(file)
            
    #%% Get files of interets
    files_of_interest = list()
    task = 'SEF' # change depending on task of interest
    matches = [task, 'proc-preproc']
    for name in fif_files:
        if all([x in name for x in matches]):
            files_of_interest.append(name)
            
    #files_of_interest.sort(key=lambda f: int(re.sub('\D','',f)))
    
    #%% Check if subject has been previously processed
    data_output = subjects_dir + '/' + task + '_MEG_output_KT'
    ecd_output = subjects_dir + '/' + task + '_ECD_output_KT'
    
    #if(op.exists(ecd_output)):
    #    continue
    
    #%% Load each file of interest
    if(len(files_of_interest) == 0):
        print('No task for',Subject)
    else:
        for file in range(0,len(files_of_interest),1):
            sample_data = Patient_files + '/' + files_of_interest[file]
            print(files_of_interest[file])
            
            #%% Check if file is in dictionary
            #chk_file = 0
            for key in my_dict:
                if key == sample_data:
                    #chk_file = 1
                    #print('key in dictionary')
                    
                    #%% Make output directories
                    if not(op.exists(data_output)):
                        os.mkdir(data_output)
                        
                    if not (op.exists(ecd_output)):
                        os.mkdir(ecd_output)
                        
                    #%% Load raw data
                    meg_preprocessed = mne.io.read_raw_fif(sample_data,preload=True)
                    
                    #%% Pick channels of interest
                    meg_preprocessed = meg_preprocessed.pick_types(meg=True,eeg=False,stim=True,eog=True,ecg=True)

                    #%% String split
                    file_details = files_of_interest[file].split('_',4)
                    
                    #%% Filter and set variables
                    if 'SEF' in file_details[2]:
                        meg_preprocessed.filter(l_freq=2,h_freq=120,phase='zero-double')
                        tmin = -0.05
                        tmax = 0.25 # 0.25 for full response, but want to exclude late field activity from dipole analysis
                        cov_max = 0.
                        baseline = (None, 0.)
                    elif 'MEF' in file_details[2]:
                        meg_preprocessed.filter(l_freq=0.1,h_freq=40,phase='zero-double')
                        tmin = -1.0
                        tmax = 0.5
                        cov_max = -0.68
                        baseline = (None, -0.68)
                    elif 'lan' in file_details[2]:
                        meg_preprocessed.filter(l_freq=1,h_freq=40,phase='zero-double')
                        tmin = -0.2
                        tmax = 0.8
                        cov_max = 0.
                        baseline = (None, 0.)
                        
                    #%% Make epochs
                    stim_channel = my_dict[key]['Stim_Channel']
                    trigger = my_dict[key]['Trigger_Value']
                    event_dict = {file_details[2]:trigger}
                    
                    events = mne.find_events(meg_preprocessed,stim_channel=stim_channel,consecutive=False)
                    meg_epochs = mne.Epochs(meg_preprocessed,events,event_id=event_dict,tmin=tmin,tmax=tmax,baseline=baseline,reject=None,
                                        preload=True,reject_by_annotation=True)
                    
                    #%% Generate and save evoked responses
                    meg_evoked = meg_epochs.average()
                    meg_evoked_topo = mne.viz.plot_evoked_topo(meg_evoked)
                    meg_evoked_plot = meg_evoked.plot(picks='meg',spatial_colors=True,gfp=True)
                    meg_evoked_plot_file = (data_output + '/' + file_details[0] + '_' + file_details[1] + '_' + file_details[2] + '_' 
                                       + file_details[3] + '_meg-evoked-plot.svg')
                    meg_evoked_plot.savefig(meg_evoked_plot_file,format='svg',transparent=True,dpi=my_dpi*10)
                    plt.close('all')
                    
                    #%% Calculate RMS of primary response
                    x = meg_evoked.copy().crop(tmin=0.015,tmax=0.060).pick_types(meg='grad') # return for crop time
                    test = x._data
                    squared_test = np.square(test)
                    mean_squared_test = np.mean(squared_test,axis=0)
                    rms = np.sqrt(mean_squared_test)
                    rms = np.reshape(rms,newshape=(1,len(rms)),order='C')
                    
                    #%% Find time of max RMS
                    times = x._raw_times
                    times = np.reshape(times,newshape=(1,len(times)),order='C')
                    evoked_max = np.argmax(rms,axis=1)
                    evoked_max_time = times[0,evoked_max]
                    evoked_max_time = np.reshape(evoked_max_time,newshape=(1,1),order='C')
                    evoked_max_time = evoked_max_time[0,0]
                    
                    #%% Compute covariance
                    meg_rank_temp = mne.compute_rank(meg_preprocessed,rank='info',proj=True)
                    meg_rank_temp = meg_rank_temp['meg']
                    rank_reduction = my_dict[key]['ICA_Removed']
                    rank = meg_rank_temp - rank_reduction
                    meg_rank = {'meg':rank}
                    meg_noise_cov = mne.compute_covariance(meg_epochs,method='auto',tmin=tmin,tmax=cov_max,rank=meg_rank)
                    
                    #%% Get transformation file
                    trans_file = my_dict[key]['Transform']
                    
                    #%% Compute surface based source space
                    src = mne.setup_source_space(subject,spacing='ico5',surface='pial',subjects_dir=subjects_dir,add_dist=True,n_jobs=-1)
                    
                    #%% Compute MEG BEM
                    meg_conductivity = (0.3,)
                    meg_model = mne.make_bem_model(subject='sample',ico=None,conductivity=meg_conductivity,subjects_dir=subjects_dir)
                    meg_bem = mne.make_bem_solution(meg_model)
                    
                    #%% Compute MEG forward operator
                    meg_fwd = mne.make_forward_solution(meg_evoked.info,trans_file,src,meg_bem,meg=True,eeg=False,mindist=5.0,n_jobs=-1)
                    meg_fwd = mne.convert_forward_solution(meg_fwd,surf_ori=True,force_fixed=False)
                    
                    #%% Compute MEG inverse operator
                    meg_inv = mne.minimum_norm.make_inverse_operator(meg_evoked.info,meg_fwd,meg_noise_cov,loose=0.2,depth=0.8,fixed=False,rank=meg_rank)
                    
                    #%% Estimate MEG SNR
                    meg_snr_white,meg_snr_est = mne.minimum_norm.estimate_snr(meg_evoked,meg_inv)
                    
                    #%% Apply inverse operator
                    method = 'sLORETA'
                    snr = 3.
                    lambda2 = 1./snr**2
                    
                    meg_stc,meg_variance = mne.minimum_norm.apply_inverse(meg_evoked,meg_inv,lambda2,method)
                    
                    #%% P20 dipole
                    evoked_ecd_P20 = meg_evoked.copy().crop(0.015,0.025)
                    dip_20 = mne.fit_dipole(evoked_ecd_P20,meg_noise_cov,meg_bem,trans=trans_file,rank=meg_rank)[0]
                    dip_20_3D = dip_20.plot_locations(trans_file,'sample',subjects_dir,mode='orthoview',coord_frame='mri',show_all=False)
                    
                    color = ['k']*len(dip_20)
                    color[np.argmax(dip_20.gof)] = 'r'
                    dip_20_Outlines = dip_20.plot_locations(trans_file,'sample',subjects_dir,mode='outlines',color=color)
                    
                    mri_20_pos = dip_20.to_mri(subject=subject,trans=trans_file,subjects_dir=subjects_dir)
                    
                    # find best fitted dipole
                    best_dip_20_idx = dip_20.gof.argmax()
                    dip_20_pos = mri_20_pos[best_dip_20_idx]
                    
                    # Grab GOF for best dipole
                    temp_20 = dip_20.gof
                    dip_20_gof = temp_20[best_dip_20_idx]
                    
                    dip_20_3D_file = (ecd_output + '/' + file_details[0] + '_' + file_details[1] + '_' + file_details[2] + '_' 
                                      + file_details[3] + '_dip-20-3D.svg')
                    dip_20_Outlines_file = (ecd_output + '/' + file_details[0] + '_' + file_details[1] + '_' + file_details[2] + '_'
                                            + file_details[3] + '_dip-20-Outlines.svg')
                    dip_20_pos_file = (ecd_output + '/' + file_details[0] + '_' + file_details[1] + '_' + file_details[2] + '_' 
                                       + file_details[3] + '_dip-20-pos.mat')
                    dip_20_gof_file = (ecd_output + '/' + file_details[0] + '_' + file_details[1] + '_' + file_details[2] + '_' 
                                       + file_details[3] + '_dip-20-gof.mat')
                    
                    dip_20_pos = np.reshape(dip_20_pos,newshape=(1,3),order='C')
                    
                    dip_20_3D.savefig(dip_20_3D_file,format='svg',transparent=True,dpi=my_dpi*10)
                    dip_20_Outlines.savefig(dip_20_Outlines_file,format='svg',transparent=True,dpi=my_dpi*10)
                    savemat(dip_20_pos_file,{'Dip_20_Pos':dip_20_pos})
                    savemat(dip_20_gof_file,{'Dip_20_GOF':dip_20_gof})
                    
                    plt.close('all')
                    
                    #%% P40 dipole
                    evoked_ecd_P40 = meg_evoked.copy().crop(0.035,0.045)
                    dip_40 = mne.fit_dipole(evoked_ecd_P40,meg_noise_cov,meg_bem,trans=trans_file,rank=meg_rank)[0]
                    dip_40_3D = dip_40.plot_locations(trans_file,'sample',subjects_dir,mode='orthoview',coord_frame='mri',show_all=False)
                    
                    color = ['k']*len(dip_40)
                    color[np.argmax(dip_40.gof)] = 'r'
                    dip_40_Outlines = dip_40.plot_locations(trans_file,'sample',subjects_dir,mode='outlines',color=color)
                    
                    mri_40_pos = dip_40.to_mri(subject=subject,trans=trans_file,subjects_dir=subjects_dir)
                    
                    # find best fitted dipole
                    best_dip_40_idx = dip_40.gof.argmax()
                    dip_40_pos = mri_40_pos[best_dip_40_idx]
                    
                    # Grab GOF for best dipole
                    temp_40 = dip_40.gof
                    dip_40_gof = temp_40[best_dip_40_idx]
                    
                    dip_40_3D_file = (ecd_output + '/' + file_details[0] + '_' + file_details[1] + '_' + file_details[2] + '_' 
                                      + file_details[3] + '_dip-40-3D.svg')
                    dip_40_Outlines_file = (ecd_output + '/' + file_details[0] + '_' + file_details[1] + '_' + file_details[2] + '_'
                                            + file_details[3] + '_dip-40-Outlines.svg')
                    dip_40_pos_file = (ecd_output + '/' + file_details[0] + '_' + file_details[1] + '_' + file_details[2] + '_' 
                                       + file_details[3] + '_dip-40-pos.mat')
                    dip_40_gof_file = (ecd_output + '/' + file_details[0] + '_' + file_details[1] + '_' + file_details[2] + '_' 
                                       + file_details[3] + '_dip-40-gof.mat')
                    
                    dip_40_pos = np.reshape(dip_40_pos,newshape=(1,3),order='C')
                    
                    dip_40_3D.savefig(dip_40_3D_file,format='svg',transparent=True,dpi=my_dpi*10)
                    dip_40_Outlines.savefig(dip_40_Outlines_file,format='svg',transparent=True,dpi=my_dpi*10)
                    savemat(dip_40_pos_file,{'Dip_40_Pos':dip_40_pos})
                    savemat(dip_40_gof_file,{'Dip_40_GOF':dip_40_gof})
                    
                    plt.close('all')
                    
                    #%% PMax dipole
                    evoked_ecd_max = meg_evoked.copy().crop(evoked_max_time,evoked_max_time)
                    dip_max = mne.fit_dipole(evoked_ecd_max,meg_noise_cov,meg_bem,trans=trans_file,rank=meg_rank)[0]
                    dip_max_3D = dip_max.plot_locations(trans_file,'sample',subjects_dir,mode='orthoview',coord_frame='mri')
                    
                    color = ['k']*len(dip_max)
                    color[np.argmax(dip_max.gof)] = 'r'
                    dip_max_Outlines = dip_max.plot_locations(trans_file,'sample',subjects_dir,mode='outlines',color=color)
                    
                    mri_max_pos = dip_max.to_mri(subject=subject,trans=trans_file,subjects_dir=subjects_dir)
                    
                    # find best fitted dipole
                    best_dip_max_idx = dip_max.gof.argmax()
                    dip_max_pos = mri_max_pos[best_dip_max_idx]
                    
                    # Grab GOF for best dipole
                    temp_max = dip_max.gof
                    dip_max_gof = temp_max[best_dip_max_idx]
                    
                    dip_max_3D_file = (ecd_output + '/' + file_details[0] + '_' + file_details[1] + '_' + file_details[2] + '_' 
                                      + file_details[3] + '_dip-max-3D.svg')
                    dip_max_Outlines_file = (ecd_output + '/' + file_details[0] + '_' + file_details[1] + '_' + file_details[2] + '_'
                                            + file_details[3] + '_dip-max-Outlines.svg')
                    dip_max_pos_file = (ecd_output + '/' + file_details[0] + '_' + file_details[1] + '_' + file_details[2] + '_' 
                                       + file_details[3] + '_dip-max-pos.mat')
                    dip_max_gof_file = (ecd_output + '/' + file_details[0] + '_' + file_details[1] + '_' + file_details[2] + '_' 
                                       + file_details[3] + '_dip-max-gof.mat')
                    
                    dip_max_pos = np.reshape(dip_max_pos,newshape=(1,3),order='C')
                    
                    dip_max_3D.savefig(dip_max_3D_file,format='svg',transparent=True,dpi=my_dpi*10)
                    dip_max_Outlines.savefig(dip_max_Outlines_file,format='svg',transparent=True,dpi=my_dpi*10)
                    savemat(dip_max_pos_file,{'Dip_Max_Pos':dip_max_pos})
                    savemat(dip_max_gof_file,{'Dip_Max_GOF':dip_max_gof})
                    
                    plt.close('all')
                    
                    #%% Extract data
                    # labels include Destrieux ('aparc.a2009s'), Desikan-Killany ('aparc'), and DKT ('DKTatlas40')
                    labels_lh = mne.read_labels_from_annot('sample',parc='aparc.a2009s',subjects_dir=subjects_dir,hemi='lh',surf_name='pial')
                    labels_rh = mne.read_labels_from_annot('sample',parc='aparc.a2009s',subjects_dir=subjects_dir,hemi='rh',surf_name='pial')
                    meg_label_lh_ts = mne.extract_label_time_course(meg_stc,labels_lh,src,mode='mean',return_generator=False)
                    meg_label_rh_ts = mne.extract_label_time_course(meg_stc,labels_rh,src,mode='mean',return_generator=False)
                    label_lh_names = [label.name for label in labels_lh]
                    label_rh_names = [label.name for label in labels_rh]
                    
                    #%% Get MEG STC locations per hemisphere and their locations
                    meg_lh_verts = meg_stc.lh_vertno
                    meg_lh_verts = meg_lh_verts.tolist()

                    # Find brain location for each vertex in the left hemisphere
                    meg_lh_verts_location = []
                    for i in range(len(meg_lh_verts)):
                        for j in range(len(labels_lh)):
                            if meg_lh_verts[i] in labels_lh[j].vertices:
                                meg_lh_verts_location.append(labels_lh[j].name)
                                
                    meg_rh_verts = meg_stc.rh_vertno
                    meg_rh_verts = meg_rh_verts.tolist()
                                
                    # Find brain location for each vertex in the right hemisphere
                    meg_rh_verts_location = []
                    for k in range(len(meg_rh_verts)):
                        for l in range(len(labels_rh)):
                            if meg_rh_verts[k] in labels_rh[l].vertices:
                                meg_rh_verts_location.append(labels_rh[l].name)
                                
                    #%% MEG source data
                    meg_lh_verts_data = meg_stc.lh_data
                    meg_rh_verts_data = meg_stc.rh_data
                    
                    #%% Convert vertex of MEG STC to MNI305 (in mm)
                    meg_lh_verts = meg_stc.lh_vertno # re-assign as was converted to list previously
                    meg_left_vertex = mne.vertex_to_mni(meg_lh_verts,0,'sample',subjects_dir=subjects_dir)
                    
                    meg_rh_verts = meg_stc.rh_vertno
                    meg_right_vertex = mne.vertex_to_mni(meg_rh_verts,1,'sample',subjects_dir=subjects_dir)
                    
                    #%% Get MEG channel names and positions
                    meg_temp = meg_preprocessed.copy().pick_types(meg=True)
                    meg_ch_pos = meg_temp._get_channel_positions(picks='meg')
                    meg_ch_names = meg_temp.ch_names
                    
                    #%% Get evoked data
                    
                    #%% Save MEG files
                    meg_epochs_file = (data_output + '/' + file_details[0] + '_' + file_details[1] + '_' + file_details[2] + '_' 
                                       + file_details[3] + '_epo.fif')
                    meg_evoked_file = (data_output + '/' + file_details[0] + '_' + file_details[1] + '_' + file_details[2] + '_'
                                       + file_details[3] + '_ave.fif')
                    meg_covariance_file = (data_output + '/' + file_details[0] + '_' + file_details[1] + '_' + file_details[2] + '_'
                                           + file_details[3] + '_cov.fif')
                    meg_source_space_file = (data_output + '/' + file_details[0] + '_' + file_details[1] + '_' + file_details[2] + '_'
                                             + file_details[3] + '_src.fif')
                    meg_bem_file = (data_output + '/' + file_details[0] + '_' + file_details[1] + '_' + file_details[2] + '_'
                                    + file_details[3] + '_bem-sol.fif')
                    meg_fwd_file = (data_output + '/' + file_details[0] + '_' + file_details[1] + '_' + file_details[2] + '_'
                                    + file_details[3] + '_fwd.fif') # NOTE: fwd solution will needs to convert to fixed orientation when loaded
                    meg_inv_file = (data_output + '/' + file_details[0] + '_' + file_details[1] + '_' + file_details[2] + '_'
                                    + file_details[3] + '_inv.fif')
                    meg_stc_file = (data_output + '/' + file_details[0] + '_' + file_details[1] + '_' + file_details[2] + '_'
                                    + file_details[3] + '_stc.stc')
                    
                    meg_epochs.save(meg_epochs_file,fmt='double',overwrite=True) # save epochs
                    meg_evoked.save(meg_evoked_file,overwrite=True) # save evoked data at sensor level
                    mne.write_cov(meg_covariance_file,meg_noise_cov,overwrite=True) # save covariance
                    mne.write_source_spaces(meg_source_space_file,src,overwrite=True) # save source space
                    mne.write_bem_solution(meg_bem_file,meg_bem,overwrite=True) # save bem solution
                    mne.write_forward_solution(meg_fwd_file,meg_fwd,overwrite=True) # save fwd solution
                    mne.minimum_norm.write_inverse_operator(meg_inv_file,meg_inv,overwrite=True) # save inv solution
                    meg_stc.save(meg_stc_file,ftype='stc',overwrite=True) # save stc
                    
                    meg_variance_file = (data_output + '/' + file_details[0] + '_' + file_details[1] + '_' + file_details[2] + '_'
                                         + file_details[3] + '_variance-explained.mat')
                    meg_snr_file = (data_output + '/' + file_details[0] + '_' + file_details[1] + '_' + file_details[2] + '_'
                                    + file_details[3] + '_snr.mat')
                    meg_evoked_position_file = (data_output + '/' + file_details[0] + '_' + file_details[1] + '_' + file_details[2] + '_'
                                                + file_details[3] + '_evoked-position.mat')
                    meg_lh_verts_data_file = (data_output + '/' + file_details[0] + '_' + file_details[1] + '_' + file_details[2] + '_'
                                              + file_details[3] + '_lh-source-mne-data.mat')
                    meg_rh_verts_data_file = (data_output + '/' + file_details[0] + '_' + file_details[1] + '_' + file_details[2] + '_'
                                              + file_details[3] + '_rh-source-mne-data.mat')
                    meg_lh_verts_position_file = (data_output + '/' + file_details[0] + '_' + file_details[1] + '_' + file_details[2] + '_'
                                                  + file_details[3] + '_lh-source-mne-position.mat')
                    meg_rh_verts_position_file = (data_output + '/' + file_details[0] + '_' + file_details[1] + '_' + file_details[2] + '_'
                                                  + file_details[3] + '_rh-source-mne-position.mat')
                    meg_lh_verts_location_file = (data_output + '/' + file_details[0] + '_' + file_details[1] + '_' + file_details[2] + '_'
                                                  + file_details[3] + '_lh-verts-location.txt')
                    meg_rh_verts_location_file = (data_output + '/' + file_details[0] + '_' + file_details[1] + '_' + file_details[2] + '_'
                                                  + file_details[3] + '_rh-verts-location.txt')
                    meg_ch_names_file = (data_output + '/' + file_details[0] + '_' + file_details[1] + '_' + file_details[2] + '_'
                                         + file_details[3] + '_ch-names.txt')
                    
                    savemat(meg_variance_file,{'Variance':meg_variance}) # variance explained by the inverse solution
                    savemat(meg_snr_file,{'SNR':meg_snr_white}) # snr estimated from the inverse solution
                    savemat(meg_evoked_position_file,{'Sensor_Position':meg_ch_pos}) # MEG channel positions for evoked data
                    #savemat(meg_source_space_file,{'Brain_Region_TS':meg_label_ts}) # average brain region activation
                    savemat(meg_lh_verts_data_file,{'LH_Verts_Data':meg_lh_verts_data}) # STC data at each source location in the left hemisphere
                    savemat(meg_rh_verts_data_file,{'RH_Verts_Data':meg_rh_verts_data}) # STC data at each source location in the right hemisphere
                    savemat(meg_lh_verts_position_file,{'LH_Verts_Position':meg_left_vertex}) # MNI305 x,y,z coordinates for each source location in the left hemisphere
                    savemat(meg_rh_verts_position_file,{'RH_Verts_Position':meg_right_vertex}) # MNI305 x,y,z coordinates for each source location in the right hemisphere
                    
                    with open(meg_ch_names_file,'w') as output:
                        for row in meg_ch_names:
                            s = ''.join(map(str,row))
                            output.write(s+'\n')
                            
                    with open(meg_lh_verts_location_file,'w') as output:
                        for row in meg_lh_verts_location:
                            s = ''.join(map(str,row))
                            output.write(s+'\n')
                            
                    with open(meg_rh_verts_location_file,'w') as output:
                        for row in meg_rh_verts_location:
                            s = ''.join(map(str,row))
                            output.write(s+'\n')
                            
                    #%% Delete variables 
                    del dip_20,dip_20_3D,dip_20_3D_file,dip_20_Outlines,dip_20_Outlines_file,dip_20_pos,dip_20_pos_file
                    del dip_40,dip_40_3D,dip_40_3D_file,dip_40_Outlines,dip_40_Outlines_file,dip_40_pos,dip_40_pos_file
                    del dip_max,dip_max_3D,dip_max_3D_file,dip_max_Outlines,dip_max_Outlines_file,dip_max_pos,dip_max_pos_file
                    del events,evoked_ecd_max,evoked_ecd_P20,evoked_ecd_P40,evoked_max,evoked_max_time,i,j,test
                    del label_lh_names,meg_bem,meg_bem_file,meg_ch_names,meg_ch_names_file,meg_ch_pos,meg_covariance_file
                    del meg_epochs,meg_epochs_file,meg_evoked,meg_evoked_file,meg_evoked_plot,meg_evoked_plot_file
                    del meg_evoked_position_file,meg_evoked_topo,meg_fwd,meg_fwd_file,meg_inv,meg_inv_file
                    del meg_left_vertex,meg_lh_verts,meg_lh_verts_data,meg_lh_verts_data_file,meg_lh_verts_location
                    del meg_lh_verts_location_file,meg_lh_verts_position_file,meg_model,meg_noise_cov,meg_preprocessed,meg_rank
                    del meg_rank_temp,meg_rh_verts,meg_rh_verts_data,meg_rh_verts_data_file,meg_rh_verts_location
                    del meg_rh_verts_location_file,meg_rh_verts_position_file,meg_right_vertex,meg_snr_est,meg_snr_file
                    del meg_snr_white,meg_source_space_file,meg_stc,meg_stc_file,meg_temp,meg_variance,meg_variance_file
                    del mri_20_pos,mri_40_pos,mri_max_pos,output,rank,rank_reduction,row,src,stim_channel
                    del trans_file,x,best_dip_20_idx,best_dip_40_idx,best_dip_max_idx,k,l,label_rh_names
                    del temp_20,temp_40,temp_max,dip_20_gof,dip_40_gof,dip_max_gof,dip_20_gof_file,dip_40_gof_file,dip_max_gof_file
                    
                    #%% Garbage collector
                    gc.collect()
                            
            #%% If file is not in dictionary (continuation from above)
            #if chk_file == 0:
            #    print(Subject,'not found in dictionary')
                
#%% Finish
print('All Done!')
