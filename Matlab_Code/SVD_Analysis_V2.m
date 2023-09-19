function output = SVD_Analysis_V2(lh_data,rh_data,lh_pos,rh_pos,lh_location,rh_location,p20,p40,pmax,p20_gof,p40_gof,pmax_gof,sub_string,side)

%% Function output = SVD_Analysis(lh_data,rh_data,lh_pos,rh_pos,sub_string)
% This function takes the path of four files from the MNE output for the
% analysis of MEG data; lh data (source localized), rh data (source
% localized), lh pos (vertex x, y, z postions), and rh pos (vertex x, y, z
% positions), as well as the location of the P20m, P40m, and Pmax dipoles, 
% and performs the following analysis: SVD analysis, percentage
% of correct points (with random permutation), percent overlap between
% dipole and distributed source (with random permutation), creates figures,
% and determines distances between dipoles and distributed sources.

%% Double check to make sure all files exist
TF = isfile(lh_data) && isfile(rh_data) && isfile(lh_pos) && isfile(rh_pos) && isfile(lh_location) && isfile(rh_location) &&...
    isfile(p20) && isfile(p40) && isfile(p20_gof) && isfile(p40_gof) && isfile(pmax_gof) && isfile(pmax);

%% Determine if subject is over the age of 19
sub_string = convertCharsToStrings(sub_string);
excluded_subjects = ["sub-unmc1006","sub-unmc1030","sub-unmc1034","sub-unmc1035","sub-unmc1108","sub-unmc1109","sub-unmc1115","sub-unmc1117",...
    "sub-unmc1127","sub-unmc1154","sub-unmc1174","sub-unmc1209","sub-unmc1229","sub-unmc1264"];

TF1 = strcmp(excluded_subjects,sub_string);
if any(TF1 == 1)
    output = 0;
    return
end

%% If all files exist, perform the analysis
if (TF == 1)

    %% Load data
    LH_Verts_Data = load(lh_data);
    RH_Verts_Data = load(rh_data);
    LH_Verts_Position = load(lh_pos);
    RH_Verts_Position = load(rh_pos);
    LH_Verts_Location = readcell(lh_location,'Delimiter',':');
    LH_Verts_Location = regexprep(LH_Verts_Location,'_','-');
    RH_Verts_Location = readcell(rh_location,'Delimiter',':');
    RH_Verts_Location = regexprep(RH_Verts_Location,'_','-');
    Dip_20_Pos = load(p20);
    Dip_40_Pos = load(p40);
    Dip_Max_Pos = load(pmax);
    Dip_20_GOF = load(p20_gof);
    Dip_40_GOF = load(p40_gof);
    Dip_Max_GOF = load(pmax_gof);

    %% Grab data from structures
    LH_Verts_Data = LH_Verts_Data.LH_Verts_Data;
    RH_Verts_Data = RH_Verts_Data.RH_Verts_Data;
    LH_Verts_Position = LH_Verts_Position.LH_Verts_Position;
    RH_Verts_Position = RH_Verts_Position.RH_Verts_Position;
    Dip_20_Pos = Dip_20_Pos.Dip_20_Pos;
    Dip_40_Pos = Dip_40_Pos.Dip_40_Pos;
    Dip_Max_Pos = Dip_Max_Pos.Dip_Max_Pos;
    Dip_20_GOF = Dip_20_GOF.Dip_20_GOF;
    Dip_40_GOF = Dip_40_GOF.Dip_40_GOF;
    Dip_Max_GOF = Dip_Max_GOF.Dip_Max_GOF;

    %% Concatenate Data
    Verts_Data = [LH_Verts_Data;RH_Verts_Data];
    Verts_Position = [LH_Verts_Position;RH_Verts_Position];
    Verts_Location = [LH_Verts_Location;RH_Verts_Location];

    %% Grab data for SVD
    SVD_Data = Verts_Data(:,66:111);

    %% Perform SVD
    [U,~,~] = svd(SVD_Data,'econ');
    % V = V';
    sef = abs(U(:,1));

    %% Find max V values
    % [val,time] = max(abs(V(1,:)));

    %% Perform spatial filtering
    filter_result = spatial_filter(sef,Verts_Position,5);

    %% Threshold filter results
    I = filter_result > 0.8*max(filter_result);
    % filter_result(~I) = 0;

    %% Grab positions of SVD identified vertices
    svd_pos = Verts_Position(I,:);

    %% Find regions associated with threshold points
    regions = Verts_Location(I);

    %% Calculate percentage of anatomically correct points
    J = contains(regions,'G-postcentral') | contains(regions,'S-central') | contains(regions,'G-precentral');
    algo_correct = (sum(J)/length(J))*100; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Find all vertex locations within 10 mm of dipole
    P20_Idx = rangesearch(Verts_Position,Dip_20_Pos,10);
    P40_Idx = rangesearch(Verts_Position,Dip_40_Pos,10);
    PMax_Idx = rangesearch(Verts_Position,Dip_Max_Pos,10);

    %% Grab source locations within 10 mm of dipole
    P20_Loc = Verts_Location(P20_Idx{1,1});
    P40_Loc = Verts_Location(P40_Idx{1,1});
    PMax_Loc = Verts_Location(PMax_Idx{1,1});

    %% Determine if distance is less than 10 mm for appropriate side
    if contains(side,'ul')
        A = contains(P20_Loc,'S-central-rh') | contains(P20_Loc,'G-precentral-rh') | contains(P20_Loc,'G-postcentral-rh');
        B = contains(P40_Loc,'S-central-rh') | contains(P40_Loc,'G-precentral-rh') | contains(P40_Loc,'G-postcentral-rh');
        C = contains(PMax_Loc,'S-central-rh') | contains(PMax_Loc,'G-precentral-rh') | contains(PMax_Loc,'G-postcentral-rh');
    elseif contains(side,'ur')
        A = contains(P20_Loc,'S-central-lh') | contains(P20_Loc,'G-precentral-lh') | contains(P20_Loc,'G-postcentral-lh');
        B = contains(P40_Loc,'S-central-lh') | contains(P40_Loc,'G-precentral-lh') | contains(P40_Loc,'G-postcentral-lh');
        C = contains(PMax_Loc,'S-central-lh') | contains(PMax_Loc,'G-precentral-lh') | contains(PMax_Loc,'G-postcentral-lh');
    end

    %% Grab vertex positions within 10 mm of dipoles
    P20_Sphere = Verts_Position(P20_Idx{1,1},:);
    P40_Sphere = Verts_Position(P40_Idx{1,1},:);
    PMax_Sphere = Verts_Position(PMax_Idx{1,1},:);

    %% Find overlap points between SVD points and dipole volumes
    x = ismember(P20_Sphere,svd_pos,'rows');
    y = ismember(P40_Sphere,svd_pos,'rows');
    z = ismember(PMax_Sphere,svd_pos,'rows');

    P20_Overlap = (sum(x)/length(x))*100; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P40_Overlap = (sum(y)/length(y))*100; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PMax_Overlap = (sum(z)/length(z))*100; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Calculate distance between each dipole and the distributed source
    P20_Dist = calculate_distance(Dip_20_Pos,svd_pos);
    P40_Dist = calculate_distance(Dip_40_Pos,svd_pos);
    PMax_Dist = calculate_distance(Dip_Max_Pos,svd_pos);

    P20_Dist = mean(P20_Dist); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P40_Dist = mean(P40_Dist); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PMax_Dist = mean(PMax_Dist); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Randomly permute area over the cortical surface
    n_iter = 10000;
    random_perm_correct = zeros(n_iter,1);
    Random_P20_Overlap = zeros(n_iter,1);
    Random_P40_Overlap = zeros(n_iter,1);
    Random_PMax_Overlap = zeros(n_iter,1);
    for random = 1:n_iter

        %% Find randomly correct points
        ix = randperm(length(Verts_Data));
        ix = ix(1:length(regions));
        random_location = Verts_Location(ix);
        random_position = Verts_Position(ix,:);
        random_correct = contains(random_location,'G-postcentral') | contains(random_location,'S-central') |...
            contains(random_location,'G-precentral');
        random_perm_correct(random,1) = (sum(random_correct)/length(random_correct))*100; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% Find random overlap with dipole volumes
        random_p20_overlap = ismember(random_position,P20_Sphere,'rows');
        random_p40_overlap = ismember(random_position,P40_Sphere,'rows');
        random_pmax_overlap = ismember(random_position,PMax_Sphere,'rows');

        Random_P20_Overlap(random,1) = (sum(random_p20_overlap)/length(random_p20_overlap))*100; %%%%%%%%%%%%%%%%%%%%%%%%%%%
        Random_P40_Overlap(random,1) = (sum(random_p40_overlap)/length(random_p40_overlap))*100; %%%%%%%%%%%%%%%%%%%%%%%%%%%
        Random_PMax_Overlap(random,1) = (sum(random_pmax_overlap)/length(random_pmax_overlap))*100; %%%%%%%%%%%%%%%%%%%%%%%%

        clear ix random_location random_position random_correct random_p20_overlap random_p40_overlap random_pmax_overlap

    end

    %% Prep for regression analysis
    if contains(side,'ul')
        J = contains(Verts_Location,'S-central-rh') | contains(Verts_Location,'G-precentral-rh') | ...
            contains(Verts_Location,'G-postcentral-rh');
    elseif contains(side,'ur')
        J = contains(Verts_Location,'S-central-lh') | contains(Verts_Location,'G-precentral-lh') | ...
            contains(Verts_Location,'G-postcentral-lh');
    end
    
    brain_regions = Verts_Location(J,1);
    vert_pos = Verts_Position(J,:);
    loc = contains(brain_regions,'G');

    %% Linear regression
    mdl = fitglm(vert_pos,loc,'linear','Distribution','binomial');

    %% Coefficients
    coeff = mdl.Coefficients.Estimate(2:end,1);
    coeff = coeff/norm(coeff); % normalize coefficients

    %% Rotate about the x-axis
    M1 = [coeff(2,1),-coeff(3,1);coeff(3,1),coeff(2,1)];
    y = [0;1];
    x1 = M1\y;
    x1 = x1/norm(x1);
    Rx = [1,0,0;0,x1(1,1),-x1(2,1);0,x1(2,1),x1(1,1)];
    b = Rx*coeff;
    
    %% Rotate about the y-axis
    M2 = [b(1,1),-b(3,1);b(3,1),b(1,1)];
    x2 = M2\y;
    x2 = x2/norm(x2);
    Ry = [x2(1,1),0,-x2(2,1);0,1,0;x2(2,1),0,x2(1,1)];

    %% Calculate rotation
    R = (Ry*Rx);

    %% Transform source locations
    Verts_Position_trans = Verts_Position*R';

    %% Get X, Y, Z distance for each dipole
    Dip_Positions = [Dip_20_Pos; Dip_40_Pos; Dip_Max_Pos];
    svd_dist = NaN(3,3); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for ii = 1:3

        %% Grab dipole position
        Dip_Pos = Dip_Positions(ii,:);

        %% Transform dipole position
        Dip_Pos_trans = Dip_Pos*R';

        %% Zero dipole position
        Verts_Position_Zero = Verts_Position_trans - Dip_Pos_trans; % Dipole at (0,0,0)

        %% Grab SVD points in new space
        svd_trans_pos = Verts_Position_Zero(I,:);

        %% Calculate SVD COM
        svd_trans_com = mean(svd_trans_pos,1);
        
        %% Save SVD Dist
        svd_dist(ii,1:end) = svd_trans_com;
    end

    %% Remove dipoles if they are more than 10 mm from the anatomical region of interest
    if (sum(A)/length(A) == 0) || isnan(sum(A)/length(A))
        svd_dist(1,:) = NaN;
        P20_Overlap = NaN;
        Random_P20_Overlap = NaN;
        P20_Dist = NaN;
    end

    if (sum(B)/length(B) == 0) || isnan(sum(B)/length(B))
        svd_dist(2,:) = NaN;
        P40_Overlap = NaN;
        Random_P40_Overlap = NaN;
        P40_Dist = NaN;
    end

    if (sum(C)/length(C) == 0) || isnan(sum(C)/length(C))
        svd_dist(3,:) = NaN;
        PMax_Overlap = NaN;
        Random_PMax_Overlap = NaN;
        PMax_Dist = NaN;
    end

    %% Create path for saving data
    save_path = append('U:\shared\database\meg-ieeg-UNMC\derivatives\',sub_string,'\ses-meg01\Matlab_SEF_Data');
    if(exist(save_path,'dir') == 0)
        mkdir(save_path)
    end

    if contains(side,'ul')
        save_file = append(save_path,'\SEFul_Data.mat');
    elseif contains(side,'ur')
        save_file = append(save_path,'\SEFur_Data.mat');
    end

    %% Save data
    save(save_file,'algo_correct','P20_Overlap','P40_Overlap','PMax_Overlap','P20_Dist','P40_Dist','PMax_Dist','random_perm_correct',...
        'Random_P20_Overlap','Random_P40_Overlap','Random_PMax_Overlap','svd_dist','sub_string','side','Dip_20_GOF','Dip_40_GOF','Dip_Max_GOF','-v7.3');

    %% Send an output
    output = 1;

else
    output = 0;
    
end

end
