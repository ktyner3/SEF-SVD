%% Clear
clear;
clc;

%% Load data
data = load('C:\Users\kevin.tyner\Desktop\MATLAB\Data\2023\MEG_SEF\Final_SEF_Analysis_Data.mat');

%% Plot random overlap with anatomical areas
x = max(size(data.Random_UL_Correct),size(data.Random_UR_Correct));
random_correct = NaN(x(1,1),2);
for ii = 1:length(data.Random_UL_Correct)
    temp = data.Random_UL_Correct{ii,1};
    temp = mean(temp);
    random_correct(ii,1) = temp;
    clear temp;
end

for ii = 1:length(data.Random_UR_Correct)
    temp = data.Random_UR_Correct{ii,1};
    temp = mean(temp);
    random_correct(ii,2) = temp;
    clear temp;
end

%% Box plot of random overlap with anatomical areas
figure(1);
fig = gca;
a = boxplot(random_correct,'colors',[1 0 0;0 0 1],'widths',0.2,'positions',[ 1 1.4 ],'Symbol','');
set(a,{'linew'},{4.0});
xticks([1 1.4]);
xticklabels({'UL','UR'});
xlabel('Task','FontSize',20);
ylabel('Percent Overlap','FontSize',20);
ylim([4 9])
fig.FontSize = 40;
fig.FontName = 'Arial';
fig.FontWeight = 'bold';
title('Random Overlap with Anatomical Regions');

%% Concatenate algorithm localization results
x = max(size(data.Algo_Correct_UL),size(data.Algo_Correct_UR));
algorithm_correct = NaN(x(1,1),2);
algorithm_correct(1:length(data.Algo_Correct_UL),1) = data.Algo_Correct_UL;
algorithm_correct(1:length(data.Algo_Correct_UR),2) = data.Algo_Correct_UR;

%% Algorithm overlap with anatomical areas
figure(2);
fig = gca;
b = boxplot(algorithm_correct,'colors',[1 0 0;0 0 1],'widths',0.2,'positions',[ 1 1.4 ],'Symbol','');
set(b,{'linew'},{4.0});
xticks([1 1.4]);
xticklabels({'UL','UR'});
xlabel('Task','FontSize',20);
ylabel('Percent Overlap','FontSize',20);
fig.FontSize = 40;
fig.FontName = 'Arial';
fig.FontWeight = 'bold';
title('Algorithm Overlap with Anatomical Regions');

%% Calculate probability algorithm localizes cortical regions better than random chance
x = max(size(data.Algo_Correct_UL),size(data.Algo_Correct_UR));
algorithm_better = NaN(x(1,1),2);
for ii = 1:length(data.Algo_Correct_UL)
    if algorithm_correct(ii,1) > random_correct(ii,1)
        algorithm_better(ii,1) = 1;
    else
        algorithm_better(ii,1) = 0;
    end
end

for ii = 1:length(data.Algo_Correct_UR)
    if algorithm_correct(ii,2) > random_correct(ii,2)
        algorithm_better(ii,2) = 1;
    else
        algorithm_better(ii,2) = 0;
    end
end

x = sum(~isnan(algorithm_better(:,2)));
y = algorithm_better(:,2) == 1;
y_ul = binocdf(sum(algorithm_better(:,1)),length(algorithm_better(:,1)),0.05,'upper');
% probability algorithm localizes cortical regions better than random chance
y_ur = binocdf(sum(y),x,0.05,'upper');

%% Calculate random area overlap with dipole volume
x = max(size(data.Random_P20_Overlap_UL),size(data.Random_P20_Overlap_UR));
random_overlap = NaN(x(1,1),6);
random_dipole_overlap = cell(x(1,1),6);
random_dipole_overlap(1:length(data.Random_P20_Overlap_UL),1) = data.Random_P20_Overlap_UL;
random_dipole_overlap(1:length(data.Random_P20_Overlap_UR),2) = data.Random_P20_Overlap_UR;
random_dipole_overlap(1:length(data.Random_P40_Overlap_UL),3) = data.Random_P40_Overlap_UL;
random_dipole_overlap(1:length(data.Random_P40_Overlap_UR),4) = data.Random_P40_Overlap_UR;
random_dipole_overlap(1:length(data.Random_PMax_Overlap_UL),5) = data.Random_PMax_Overlap_UL;
random_dipole_overlap(1:length(data.Random_PMax_Overlap_UR),6) = data.Random_PMax_Overlap_UR;

for ii = 1:size(random_dipole_overlap,2)
    for jj = 1:size(random_dipole_overlap,1)
        random_overlap(jj,ii) = mean(random_dipole_overlap{jj,ii});
    end
end

%% Plot random dipole overlap
figure(3);
fig = gca;
c = boxplot(random_overlap,'colors',[1 0 0;0 0 1],'widths',0.2,'positions',[ 1 1.2 2 2.2 3 3.2],'Symbol','');
set(c,{'linew'},{4.0});
xticks([1.1 2.1 3.1]);
xlabel('Dipole','FontSize',20);
ylabel('Percent Overlap','FontSize',20);
ylim([0 1])
xticklabels({'P20', 'P40', 'PMax'});
fig.FontSize = 35;
fig.FontName = 'Arial';
fig.FontWeight = 'bold';
title('Random Location Overlap with Dipole Volume');
legend(findobj(gca,'Tag','Box'),'UR','UL')

%% Grab dipole overlap with algorithm
x = max(size(data.P20_UL_Overlap),size(data.P20_UR_Overlap));
y = max(size(data.P40_UL_Overlap),size(data.P40_UR_Overlap));
z = max(size(data.PMax_UL_Overlap),size(data.PMax_UR_Overlap));
max_size = [x(1,1),y(1,1),z(1,1)];
zz = max(max_size);
dipole_algorithm_overlap = NaN(zz,6);
dipole_algorithm_overlap(1:length(data.P20_UL_Overlap),1) = data.P20_UL_Overlap;
dipole_algorithm_overlap(1:length(data.P20_UR_Overlap),2) = data.P20_UR_Overlap;
dipole_algorithm_overlap(1:length(data.P40_UL_Overlap),3) = data.P40_UL_Overlap;
dipole_algorithm_overlap(1:length(data.P40_UR_Overlap),4) = data.P40_UR_Overlap;
dipole_algorithm_overlap(1:length(data.PMax_UL_Overlap),5) = data.PMax_UL_Overlap;
dipole_algorithm_overlap(1:length(data.PMax_UR_Overlap),6) = data.PMax_UR_Overlap;

%% Plot dipole overlap with algorithm
figure(4);
fig = gca;
d = boxplot(dipole_algorithm_overlap,'colors',[1 0 0;0 0 1],'widths',0.2,'positions',[ 1 1.2 2 2.2 3 3.2],'Symbol','');
set(d,{'linew'},{4.0});
xticks([1.1 2.1 3.1]);
xlabel('Dipole','FontSize',20);
ylabel('Percent Overlap','FontSize',20);
ylim([0 100])
xticklabels({'P20m', 'P40m', 'PMax'});
fig.FontSize = 40;
fig.FontName = 'Arial';
fig.FontWeight = 'bold';
title('Dipole Overlap with Algorithm Area');
legend(findobj(gca,'Tag','Box'),'UR','UL')

%% Calculate percentage of time random overlap is greater than algorithm overlap
percent_random_overlap_greater = NaN(zz,6);
for ii = 1:size(random_dipole_overlap,2)
    for jj = 1:size(random_dipole_overlap,1)
        temp = random_dipole_overlap{jj,ii};
        I = temp > dipole_algorithm_overlap(jj,ii);
        I = (sum(I)/length(I))*100;
        percent_random_overlap_greater(jj,ii) = I;
    end
end

%% Plot random greater
figure(5);
fig = gca;
e = boxplot(percent_random_overlap_greater,'colors',[1 0 0;0 0 1],'widths',0.2,'positions',[ 1 1.2 2 2.2 3 3.2],'Symbol','');
set(e,{'linew'},{4.0});
xticks([1.1 2.1 3.1]);
xlabel('Dipole','FontSize',20);
ylabel('Percent of Random Trials','FontSize',20);
ylim([0 25])
xticklabels({'P20m', 'P40m', 'PMax'});
fig.FontSize = 40;
fig.FontName = 'Arial';
fig.FontWeight = 'bold';
title('Random Overlap Greater than Algorithm Overlap');
legend(findobj(gca,'Tag','Box'),'UR','UL')

%% Grab dipole distance to distributed source in original space
dipole_distance = NaN(zz,6);
dipole_distance(1:length(data.P20_UL_Dist),1) = data.P20_UL_Dist;
dipole_distance(1:length(data.P20_UR_Dist),2) = data.P20_UR_Dist;
dipole_distance(1:length(data.P40_UL_Dist),3) = data.P40_UL_Dist;
dipole_distance(1:length(data.P40_UR_Dist),4) = data.P40_UR_Dist;
dipole_distance(1:length(data.PMax_UL_Dist),5) = data.PMax_UL_Dist;
dipole_distance(1:length(data.PMax_UR_Dist),6) = data.PMax_UR_Dist;

%% Plot dipole distance to algorithm localization
figure(6);
fig = gca;
f = boxplot(dipole_distance,'colors',[1 0 0;0 0 1],'widths',0.2,'positions',[ 1 1.2 2 2.2 3 3.2],'Symbol','');
set(f,{'linew'},{4.0});
xticks([1.1 2.1 3.1]);
xlabel('Dipole','FontSize',20);
ylabel('Distance to Distributed Source (mm)','FontSize',20);
ylim([0 40])
xticklabels({'P20m', 'P40m', 'PMax'});
fig.FontSize = 32;
fig.FontName = 'Arial';
fig.FontWeight = 'bold';
title('Dipole Distance from Algorithm Localization');
legend(findobj(gca,'Tag','Box'),'UR','UL')

%% Grab transformed X, Y, Z coordinates for each dipole
x = max(size(data.P20_SVD_UL),size(data.P20_SVD_UR));
y = max(size(data.P40_SVD_UL),size(data.P40_SVD_UR));
z = max(size(data.PMax_SVD_UL),size(data.PMax_SVD_UR));

P20_transformed = NaN(x(1,1),6);
P20_transformed(1:length(data.P20_SVD_UL),1:3) = data.P20_SVD_UL;
P20_transformed(1:length(data.P20_SVD_UR),4:6) = data.P20_SVD_UR;

P40_transformed = NaN(y(1,1),6);
P40_transformed(1:length(data.P40_SVD_UL),1:3) = data.P40_SVD_UL;
P40_transformed(1:length(data.P40_SVD_UR),4:6) = data.P40_SVD_UR;

PMax_transformed = NaN(z(1,1),6);
PMax_transformed(1:length(data.PMax_SVD_UL),1:3) = data.PMax_SVD_UL;
PMax_transformed(1:length(data.PMax_SVD_UR),4:6) = data.PMax_SVD_UR;

%% Make bar plots of the distribution of x, y, and z points for P20
figure(7)
t = tiledlayout(3,2);
txt = title(t,'P20');
txt.FontWeight = 'bold';
txt.FontSize = 32;
xlabel(t,'Distance (mm)','FontSize',20,'FontWeight','bold')
ylabel(t,'Number of Dipoles','FontSize',20,'FontWeight','bold')

% Tile 1
nexttile
histogram(P20_transformed(:,1),20);
title('X UL');
xlim([-30 30])
ylim([1 12])

% Tile 2
nexttile
histogram(P20_transformed(:,4),20);
title('X UR');
xlim([-30 30])
ylim([1 12])

% Tile 3
nexttile
histogram(P20_transformed(:,2),20);
title('Y UL');
xlim([-30 30])
ylim([1 12])

% Tile 4
nexttile
histogram(P20_transformed(:,5),20);
title('Y UR');
xlim([-30 30])
ylim([1 12])

% Tile 5
nexttile
histogram(P20_transformed(:,3),20);
title('Z UL');
xlim([-30 30])
ylim([1 12])

% Tile 6
nexttile
histogram(P20_transformed(:,6),20);
title('Z UR');
xlim([-30 30])
ylim([1 12])

%% Make bar plots of the distribution of x, y, and z points for P40
figure(8)
t = tiledlayout(3,2);
txt = title(t,'P40');
txt.FontWeight = 'bold';
txt.FontSize = 32;
xlabel(t,'Distance (mm)','FontSize',20,'FontWeight','bold')
ylabel(t,'Number of Dipoles','FontSize',20,'FontWeight','bold')

% Tile 1
nexttile
histogram(P40_transformed(:,1),20);
title('X UL');
xlim([-30 30])
ylim([1 12])

% Tile 2
nexttile
histogram(P40_transformed(:,4),20);
title('X UR');
xlim([-30 30])
ylim([1 12])

% Tile 3
nexttile
histogram(P40_transformed(:,2),20);
title('Y UL');
xlim([-30 30])
ylim([1 12])

% Tile 4
nexttile
histogram(P40_transformed(:,5),20);
title('Y UR');
xlim([-30 30])
ylim([1 12])

% Tile 5
nexttile
histogram(P40_transformed(:,3),20);
title('Z UL');
xlim([-30 30])
ylim([1 12])

% Tile 6
nexttile
histogram(P40_transformed(:,6),20);
title('Z UR');
xlim([-30 30])
ylim([1 12])

%% Make bar plots of the distribution of x, y, and z points for PMax
figure(9)
t = tiledlayout(3,2);
txt = title(t,'PMax');
txt.FontWeight = 'bold';
txt.FontSize = 32;
xlabel(t,'Distance (mm)','FontSize',20,'FontWeight','bold')
ylabel(t,'Number of Dipoles','FontSize',20,'FontWeight','bold')

% Tile 1
nexttile
histogram(PMax_transformed(:,1),20);
title('X UL');
xlim([-30 30])
ylim([1 12])

% Tile 2
nexttile
histogram(PMax_transformed(:,4),20);
title('X UR');
xlim([-30 30])
ylim([1 12])

% Tile 3
nexttile
histogram(PMax_transformed(:,2),20);
title('Y UL');
xlim([-30 30])
ylim([1 12])

% Tile 4
nexttile
histogram(PMax_transformed(:,5),20);
title('Y UR');
xlim([-30 30])
ylim([1 12])

% Tile 5
nexttile
histogram(PMax_transformed(:,3),20);
title('Z UL');
xlim([-30 30])
ylim([1 12])

% Tile 6
nexttile
histogram(PMax_transformed(:,6),20);
title('Z UR');
xlim([-30 30])
ylim([1 12])

%% Statistics
p_val20 = zeros(6,1);
p_val40 = zeros(6,1);
p_valMax = zeros(6,1);

for ii = 1:6
    z = P20_transformed(:,ii);
    zz = P40_transformed(:,ii);
    zzz = PMax_transformed(:,ii);

    p_val20(ii,1) = signrank(z);
    p_val40(ii,1) = signrank(zz);
    p_valMax(ii,1) = signrank(zzz);

    clear z zz zzz

end

%% Plot GOF values
x = max(size(data.Dip_20_GOF_UL),size(data.Dip_20_GOF_UR));
GOF = NaN(x(1,1),6);
GOF(1:length(data.Dip_20_GOF_UL),1) = cell2mat(data.Dip_20_GOF_UL);
GOF(1:length(data.Dip_20_GOF_UR),2) = cell2mat(data.Dip_20_GOF_UR);
GOF(1:length(data.Dip_40_GOF_UL),3) = cell2mat(data.Dip_40_GOF_UL);
GOF(1:length(data.Dip_40_GOF_UR),4) = cell2mat(data.Dip_40_GOF_UR);
GOF(1:length(data.Dip_Max_GOF_UL),5) = cell2mat(data.Dip_Max_GOF_UL);
GOF(1:length(data.Dip_Max_GOF_UR),6) = cell2mat(data.Dip_Max_GOF_UR);
figure(10)
fig = gca;
g = boxplot(GOF,'colors',[1 0 0;0 0 1],'widths',0.2,'positions',[ 1 1.2 2 2.2 3 3.2],'Symbol','');
set(g,{'linew'},{4.0});
xticks([1.1 2.1 3.1]);
xlabel('Dipole','FontSize',20);
ylabel('GOF','FontSize',20);
ylim([0 110])
xticklabels({'P20m', 'P40m', 'PMax'});
fig.FontSize = 32;
fig.FontName = 'Arial';
fig.FontWeight = 'bold';
title('GOF Values');
legend(findobj(gca,'Tag','Box'),'UR','UL')

%% Estimate false discovery rate with Benjamini-Hochberg procedure
total_stats = [p_val20;p_val40;p_valMax];
for ii = 1:length(total_stats)
    corrected(ii,1) = 0.05/ii;
end
total_stats = sort(total_stats,'descend');
overall_stats = [corrected,total_stats];

%% Calculate probability algorithm/dipole overlap is larger than chance/dipole overlap
x = max(size(data.P20_UL_Overlap),size(data.P20_UR_Overlap));
overlaps = NaN(x(1,1),6);

% P20 UL Overlap
for ii = 1:length(data.P20_UL_Overlap)
    if isnan(data.P20_UL_Overlap(ii,1))
        continue
    else
        temp = data.Random_P20_Overlap_UL{ii,1};
        temp = mean(temp);
        if data.P20_UL_Overlap(ii,1) > temp
            overlaps(ii,1) = 1;
        else
            overlaps(ii,1) = 0;
        end
    end
end

% P20 UR Overlap
for ii = 1:length(data.P20_UR_Overlap)
    if isnan(data.P20_UR_Overlap(ii,1))
        continue
    else
        temp = data.Random_P20_Overlap_UR{ii,1};
        temp = mean(temp);
        if data.P20_UR_Overlap(ii,1) > temp
            overlaps(ii,2) = 1;
        else
            overlaps(ii,2) = 0;
        end
    end
end

% P40 UL Overlap
for ii = 1:length(data.P40_UL_Overlap)
    if isnan(data.P40_UL_Overlap(ii,1))
        continue
    else
        temp = data.Random_P40_Overlap_UL{ii,1};
        temp = mean(temp);
        if data.P40_UL_Overlap(ii,1) > temp
            overlaps(ii,3) = 1;
        else
            overlaps(ii,3) = 0;
        end
    end
end

% P40 UR Overlap
for ii = 1:length(data.P40_UR_Overlap)
    if isnan(data.P40_UR_Overlap(ii,1))
        continue
    else
        temp = data.Random_P40_Overlap_UR{ii,1};
        temp = mean(temp);
        if data.P40_UR_Overlap(ii,1) > temp
            overlaps(ii,4) = 1;
        else
            overlaps(ii,4) = 0;
        end
    end
end

% PMax UL Overlap
for ii = 1:length(data.PMax_UL_Overlap)
    if isnan(data.PMax_UL_Overlap(ii,1))
        continue
    else
        temp = data.Random_PMax_Overlap_UL{ii,1};
        temp = mean(temp);
        if data.PMax_UL_Overlap(ii,1) > temp
            overlaps(ii,5) = 1;
        else
            overlaps(ii,5) = 0;
        end
    end
end

% PMax UR Overlap
for ii = 1:length(data.PMax_UR_Overlap)
    if isnan(data.PMax_UR_Overlap(ii,1))
        continue
    else
        temp = data.Random_PMax_Overlap_UR{ii,1};
        temp = mean(temp);
        if data.PMax_UR_Overlap(ii,1) > temp
            overlaps(ii,6) = 1;
        else
            overlaps(ii,6) = 0;
        end
    end
end

% Calculate probability values
p20_ul = binocdf(sum(overlaps(:,1),'omitnan'),sum(~isnan(overlaps(:,1))),0.05,'upper');
p20_ur = binocdf(sum(overlaps(:,2),'omitnan'),sum(~isnan(overlaps(:,2))),0.05,'upper');
p40_ul = binocdf(sum(overlaps(:,3),'omitnan'),sum(~isnan(overlaps(:,3))),0.05,'upper');
p40_ur = binocdf(sum(overlaps(:,4),'omitnan'),sum(~isnan(overlaps(:,4))),0.05,'upper');
pmax_ul = binocdf(sum(overlaps(:,5),'omitnan'),sum(~isnan(overlaps(:,5))),0.05,'upper');
pmax_ur = binocdf(sum(overlaps(:,6),'omitnan'),sum(~isnan(overlaps(:,6))),0.05,'upper');
