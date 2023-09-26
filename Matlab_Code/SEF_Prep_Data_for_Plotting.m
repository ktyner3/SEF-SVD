%% Clear
clear;
clc;

%% Set variables
count_ul = 0;
count_ur = 0;

%% Set directory
sub_dir = 'U:\shared\database\meg-ieeg-UNMC\derivatives\';

%% Grab participants to include in the analysis
participants = readtable('U:\shared\users\ktyner\Papers\2023\NeuroImage_Clinical\New_SEF_Demographics.xlsx');
participants = string(participants{:,1});
participants = append('sub-',participants);

%% Exclude Subjects
exclude_ul = ["sub-unmc1079","sub-unmc1152","sub-unmc1195","sub-unmc1198"];
exclude_ur = ["sub-unmc1000","sub-unmc1033","sub-unmc1061","sub-unmc1143","sub-unmc1156","sub-unmc1165","sub-unmc1192","sub-unmc1196",...
    "sub-unmc1203","sub-unmc1243","sub-unmc1247","sub-unmc1266"];

%% Loop through subjects and grab those that were analyzed
for ii = 1:length(participants)

    %% Grab subject data
    sub_file = append(sub_dir,participants(ii),'\ses-meg01\Matlab_SEF_Data\');

    %% Grab file names
    sub_data = dir(sub_file);
    sub_data(1:2,:) = [];

    %% Loop to load data
    for jj = 1:length(sub_data)

        %% Grab file name
        sub_data_file = append(sub_file,sub_data(jj).name);

        %% Load data
        file_data = load(sub_data_file);

        %% Grab data
        if contains(sub_data(jj).name,'SEFul')

            %% Check to see if subject file should be included
            TF1 = strcmp(participants(ii,1),exclude_ul);
            if any(TF1 == 1)
                continue
            else
                count_ul = count_ul + 1;
                P20_SVD_UL(count_ul,1:3) = file_data.svd_dist(1,:);
                P40_SVD_UL(count_ul,1:3) = file_data.svd_dist(2,:);
                PMax_SVD_UL(count_ul,1:3) = file_data.svd_dist(3,:);
                P20_UL_Dist(count_ul,1) = file_data.P20_Dist;
                P40_UL_Dist(count_ul,1) = file_data.P40_Dist;
                PMax_UL_Dist(count_ul,1) = file_data.PMax_Dist;
                P20_UL_Overlap(count_ul,1) = file_data.P20_Overlap;
                P40_UL_Overlap(count_ul,1) = file_data.P40_Overlap;
                PMax_UL_Overlap(count_ul,1) = file_data.PMax_Overlap;
                Algo_Correct_UL(count_ul,1) = file_data.algo_correct;
                Sub_String_UL{count_ul,1} = file_data.sub_string;
                Random_P20_Overlap_UL{count_ul,1} = file_data.Random_P20_Overlap;
                Random_P40_Overlap_UL{count_ul,1} = file_data.Random_P40_Overlap;
                Random_PMax_Overlap_UL{count_ul,1} = file_data.Random_PMax_Overlap;
                Random_UL_Correct{count_ul,1} = file_data.random_perm_correct;
                Dip_20_GOF_UL{count_ul,1} = file_data.Dip_20_GOF;
                Dip_40_GOF_UL{count_ul,1} = file_data.Dip_40_GOF;
                Dip_Max_GOF_UL{count_ul,1} = file_data.Dip_Max_GOF;
            end
            
        elseif contains(sub_data(jj).name,'SEFur')
            TF2 = strcmp(participants(ii,1),exclude_ur);
            if any(TF2 == 1)
                continue
            else
                count_ur = count_ur + 1;
                P20_SVD_UR(count_ur,1:3) = file_data.svd_dist(1,:);
                P40_SVD_UR(count_ur,1:3) = file_data.svd_dist(2,:);
                PMax_SVD_UR(count_ur,1:3) = file_data.svd_dist(3,:);
                P20_UR_Dist(count_ur,1) = file_data.P20_Dist;
                P40_UR_Dist(count_ur,1) = file_data.P40_Dist;
                PMax_UR_Dist(count_ur,1) = file_data.PMax_Dist;
                P20_UR_Overlap(count_ur,1) = file_data.P20_Overlap;
                P40_UR_Overlap(count_ur,1) = file_data.P40_Overlap;
                PMax_UR_Overlap(count_ur,1) = file_data.PMax_Overlap;
                Algo_Correct_UR(count_ur,1) = file_data.algo_correct;
                Sub_String_UR{count_ur,1} = file_data.sub_string;
                Random_P20_Overlap_UR{count_ur,1} = file_data.Random_P20_Overlap;
                Random_P40_Overlap_UR{count_ur,1} = file_data.Random_P40_Overlap;
                Random_PMax_Overlap_UR{count_ur,1} = file_data.Random_PMax_Overlap;
                Random_UR_Correct{count_ur,1} = file_data.random_perm_correct;
                Dip_20_GOF_UR{count_ur,1} = file_data.Dip_20_GOF;
                Dip_40_GOF_UR{count_ur,1} = file_data.Dip_40_GOF;
                Dip_Max_GOF_UR{count_ur,1} = file_data.Dip_Max_GOF;
            end
        end
    end
end

%% Save Data
save('C:\xxxxxxx\Final_SEF_Analysis_Data.mat','Algo_Correct_UL','Algo_Correct_UR','P20_SVD_UL',...
    'P20_SVD_UR','P20_UL_Dist','P20_UL_Overlap','P20_UR_Dist','P20_UR_Overlap','P40_SVD_UL','P40_SVD_UR','P40_UL_Dist','P40_UL_Overlap',...
    'P40_UR_Dist','P40_UR_Overlap','PMax_SVD_UL','PMax_SVD_UR','PMax_UL_Dist','PMax_UL_Overlap','PMax_UR_Dist','PMax_UR_Overlap',...
    'Random_P20_Overlap_UL','Random_P20_Overlap_UR','Random_P40_Overlap_UL','Random_P40_Overlap_UR','Random_PMax_Overlap_UL','Random_PMax_Overlap_UR',...
    'Random_UL_Correct','Random_UR_Correct','Sub_String_UL','Sub_String_UR','Dip_20_GOF_UL','Dip_40_GOF_UL','Dip_Max_GOF_UL',...
    'Dip_20_GOF_UR','Dip_40_GOF_UR','Dip_Max_GOF_UR');
