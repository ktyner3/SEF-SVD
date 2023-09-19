%% Clear
clear;
clc;

%% Add FT
%addpath('U:\shared\database\meg-ieeg-UNMC\fieldtrip\fieldtrip');
%ft_defaults;

%% Add functions
addpath('U:\shared\users\ktyner\2023\Functions\');

%% Set directory
sub_dir = 'U:\shared\database\meg-ieeg-UNMC\derivatives\';
subs = dir(sub_dir);
subs = {subs([subs.isdir]).name};
subs = subs(~ismember(subs,{'.','..'}));
subs = string(subs)';

%% Grab participants to include in the analysis
participants = readtable('U:\shared\users\ktyner\Papers\2023\NeuroImage_Clinical\New_SEF_Demographics.xlsx');
participants = string(participants{:,1});
participants = append('sub-',participants);

%% Analyze each subject
for ii = 1:length(participants)

    %% See if subject is included in the analysis
    sub_string = participants(ii,1);
    TF = strcmp(sub_string,participants);

    %% If they are, grab them
    if any(TF == 1)

        %% Go to subject folder and figure out which task(s) they have
        sub_svd_data_folder = append(sub_dir,sub_string,'\ses-meg01\SEF_MEG_output_KT\');
        sub_ecd_data_folder = append(sub_dir,sub_string,'\ses-meg01\SEF_ECD_output_KT\');

        sub_svd_files = dir(fullfile(sub_svd_data_folder,'*.mat'));
        sub_ecd_files = dir(fullfile(sub_ecd_data_folder,'*.mat'));

        svd_file_names = {sub_svd_files.name}';
        svd_file_names = string(svd_file_names);
        ecd_file_names = {sub_ecd_files.name}';
        ecd_file_names = string(ecd_file_names);

        left_svd_files = contains(svd_file_names,'SEFul');
        right_svd_files = contains(svd_file_names,'SEFur');

        left_ecd_files = contains(ecd_file_names,'SEFul');
        right_ecd_files = contains(ecd_file_names,'SEFur');

        %% Do the analysis for either task
        if any(left_svd_files == 1)
            lh_data = append(sub_svd_data_folder,sub_string,'_ses-preop_task-SEFul_run-01_lh-source-mne-data.mat');
            rh_data = append(sub_svd_data_folder,sub_string,'_ses-preop_task-SEFul_run-01_rh-source-mne-data.mat');
            lh_pos = append(sub_svd_data_folder,sub_string,'_ses-preop_task-SEFul_run-01_lh-source-mne-position.mat');
            rh_pos = append(sub_svd_data_folder,sub_string,'_ses-preop_task-SEFul_run-01_rh-source-mne-position.mat');
            lh_location = append(sub_svd_data_folder,sub_string,'_ses-preop_task-SEFul_run-01_lh-verts-location.txt');
            rh_location = append(sub_svd_data_folder,sub_string,'_ses-preop_task-SEFul_run-01_rh-verts-location.txt');
            p20 = append(sub_ecd_data_folder,sub_string,'_ses-preop_task-SEFul_run-01_dip-20-pos.mat');
            p40 = append(sub_ecd_data_folder,sub_string,'_ses-preop_task-SEFul_run-01_dip-40-pos.mat');
            pmax = append(sub_ecd_data_folder,sub_string,'_ses-preop_task-SEFul_run-01_dip-max-pos.mat');
            p20_gof = append(sub_ecd_data_folder,sub_string,'_ses-meg01_task-SEFul_run-01_dip-20-gof.mat');
            p40_gof = append(sub_ecd_data_folder,sub_string,'_ses-meg01_task-SEFul_run-01_dip-40-gof.mat');
            pmax_gof = append(sub_ecd_data_folder,sub_string,'_ses-meg01_task-SEFul_run-01_dip-max-gof.mat');

            TF1 = isfile(lh_data) && isfile(rh_data) && isfile(lh_pos) && isfile(rh_pos) && isfile(lh_location) &&...
                isfile(rh_location) && isfile(p20) && isfile(p40) && isfile(pmax) && isfile(p20_gof) && ...
                isfile(p40_gof) && isfile(pmax_gof);
            if (TF1 == 1)
                fprintf('Analyzing %s task %s..\n',sub_string,'SEFul');
                output = SVD_Analysis_V2(lh_data,rh_data,lh_pos,rh_pos,lh_location,rh_location,p20,p40,pmax,p20_gof,...
                    p40_gof,pmax_gof,sub_string,'SEFul');
            end
        end

        if any(right_svd_files == 1)
            lh_data = append(sub_svd_data_folder,sub_string,'_ses-preop_task-SEFur_run-01_lh-source-mne-data.mat');
            rh_data = append(sub_svd_data_folder,sub_string,'_ses-preop_task-SEFur_run-01_rh-source-mne-data.mat');
            lh_pos = append(sub_svd_data_folder,sub_string,'_ses-preop_task-SEFur_run-01_lh-source-mne-position.mat');
            rh_pos = append(sub_svd_data_folder,sub_string,'_ses-preop_task-SEFur_run-01_rh-source-mne-position.mat');
            lh_location = append(sub_svd_data_folder,sub_string,'_ses-preop_task-SEFul_run-01_lh-verts-location.txt');
            rh_location = append(sub_svd_data_folder,sub_string,'_ses-preop_task-SEFul_run-01_rh-verts-location.txt');
            p20 = append(sub_ecd_data_folder,sub_string,'_ses-preop_task-SEFur_run-01_dip-20-pos.mat');
            p40 = append(sub_ecd_data_folder,sub_string,'_ses-preop_task-SEFur_run-01_dip-40-pos.mat');
            pmax = append(sub_ecd_data_folder,sub_string,'_ses-preop_task-SEFur_run-01_dip-max-pos.mat');
            p20_gof = append(sub_ecd_data_folder,sub_string,'_ses-meg01_task-SEFur_run-01_dip-20-gof.mat');
            p40_gof = append(sub_ecd_data_folder,sub_string,'_ses-meg01_task-SEFur_run-01_dip-40-gof.mat');
            pmax_gof = append(sub_ecd_data_folder,sub_string,'_ses-meg01_task-SEFur_run-01_dip-max-gof.mat');

            TF1 = isfile(lh_data) && isfile(rh_data) && isfile(lh_pos) && isfile(rh_pos) && isfile(lh_location) &&...
                isfile(rh_location) && isfile(p20) && isfile(p40) && isfile(pmax) && isfile(p20_gof) && ...
                isfile(p40_gof) && isfile(pmax_gof);
            if (TF1 == 1)
                fprintf('Analyzing %s task %s..\n',sub_string,'SEFur');
                output = SVD_Analysis_V2(lh_data,rh_data,lh_pos,rh_pos,lh_location,rh_location,p20,p40,pmax,p20_gof,...
                    p40_gof,pmax_gof,sub_string,'SEFur');
            end
        end
        
    else
        continue
    end
end
