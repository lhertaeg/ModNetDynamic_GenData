%% Statistical analysis of the 'SimpAdEx' model parameters for mGluR3
% This worksheet contains all code and further explanations to perform statistical
% analyses of the model parameters for mGluR3

%% Loading and processing data
% We load the data containing all parameters

load('EverythingInOneFile34.mat');

%%
% For the subsequent analyses, we only need some of the information.
% Therefore, we generate a new data structure / data array:

pharma = reshape(pharma,[690 1]);
type = reshape(type,[690 1]);
LoV = ["genotype","age2","cell_type2","layer","species","pharma","type"];

for i=1:690
    for j=1:size(LoV,2)
        if eval("~isa(" + LoV(j) + "{" + i + ",1},'char')")
            eval(LoV(j) + "{" + i + ",1}='unknown';")
        end
    end
end

Data{:,1} = genotype;
Data{:,2} = age2;
Data{:,3} = cell_type2;
Data{:,4} = layer;
Data{:,5} = species;
Data{:,6} = pharma;
Data{:,7} = type;
Data{:,8} = C;          
Data{:,9} = gL;
Data{:,10} = tau;
Data{:,11} = EL;
Data{:,12} = sf;
Data{:,13} = Vth;
Data{:,14} = Vr;
Data{:,15} = Vup;
Data{:,16} = tcw;
Data{:,17} = b;

%% 
% We create new data frame that contains mGlu-data only

ind_WT = [ind_mGluR3WT_L3_PC',ind_mGluR3WT_L5_pp',ind_mGluR3WT_BB'];
ind_KO = [ind_mGluR3KO_L3_PC',ind_mGluR3KO_L5_nn',ind_mGluR3KO_FS'];
ind_KO_cont = [ind_mGluR3KO_L3_PC_Glu_control',ind_mGluR3KO_L5_nn_Glu_control',ind_mGluR3KO_FS_Glu_control'];
ind_KO_agon = [ind_mGluR3KO_L3_PC_Glu_agonist',ind_mGluR3KO_L5_nn_Glu_agonist',ind_mGluR3KO_FS_Glu_agonist'];
ind_WT_agon = [ind_mGluR3WT_L3_PC_Glu_agonist',ind_mGluR3WT_L5_pp_Glu_agonist',ind_mGluR3WT_BB_Glu_agonist'];
ind_WT_cont = [ind_mGluR3WT_L3_PC_Glu_control',ind_mGluR3WT_L5_pp_Glu_control',ind_mGluR3WT_BB_Glu_control'];

Data_WT = Data(:,[2:5,8:17]);
Data_WT = cellfun(@(x) x(ind_WT), Data_WT, 'UniformOutput', false);
Data_KO = Data(:,[2:5,8:17]);
Data_KO = cellfun(@(x) x(ind_KO), Data_KO, 'UniformOutput', false);
Data_WT_cont = Data(:,[2:5,8:17]);
Data_WT_cont = cellfun(@(x) x(ind_WT_cont), Data_WT_cont, 'UniformOutput', false);
Data_KO_cont = Data(:,[2:5,8:17]);
Data_KO_cont = cellfun(@(x) x(ind_KO_cont), Data_KO_cont, 'UniformOutput', false);
Data_WT_agon = Data(:,[2:5,8:17]);
Data_WT_agon = cellfun(@(x) x(ind_WT_agon), Data_WT_agon, 'UniformOutput', false);
Data_KO_agon = Data(:,[2:5,8:17]);
Data_KO_agon = cellfun(@(x) x(ind_KO_agon), Data_KO_agon, 'UniformOutput', false);
ColNam = ["age","cell_type","layer","species","C","gL","tau","EL","sf","Vth","Vr","Vup","tcw","b"];

%% Analysis of WT vs. KO data
% In a first step we look at the number of data points with respect to
% specific conditions (age, cell type, layer, species)

N_WT = size(Data_WT{1,1},1);
disp(['The WT data set contains ' num2str(N_WT) ' data points'])
class_age_WT = unique(Data_WT{1,1});
disp(['The following age categories are present for WT: ' strjoin(class_age_WT,', ')])
class_celltype_WT = unique(Data_WT{1,2});
disp(['The following cell types are present for control: ' strjoin(class_celltype_WT,', ')])
disp(['BB: ' num2str(length(ind_mGluR3WT_BB)) ', L3-PC: ' num2str(length(ind_mGluR3WT_L3_PC)) ', L5-PC: ' num2str(length(ind_mGluR3WT_L5_pp))])
class_layer_WT = unique(Data_WT{1,3});
disp(['The following layers are present for WT: ', strjoin(class_layer_WT,', ')])
class_species_WT = unique(Data_WT{1,4});
disp(['The recordings come from the follwoing species: ', strjoin(class_species_WT,', ')])
disp([' '])
N_KO = size(Data_KO{1,1},1);
disp(['The KO data set contains ' num2str(N_KO) ' data points'])
class_age_KO = unique(Data_KO{1,1});
disp(['The following age categories are present for KO: ' strjoin(class_age_KO,', ')])
class_celltype_KO = unique(Data_KO{1,2});
disp(['The following cell types are present for KO: ' strjoin(class_celltype_KO,', ')])
disp(['FS: ' num2str(length(ind_mGluR3KO_FS)) ', L3-PC: ' num2str(length(ind_mGluR3KO_L3_PC)) ', L5-PC: ' num2str(length(ind_mGluR3KO_L5_nn))])
class_layer_KO = unique(Data_KO{1,3});
disp(['The following layers are present for KO: ', strjoin(class_layer_KO,', ')])
class_species_KO = unique(Data_KO{1,4});
disp(['The recordings come from the follwoing species: ', strjoin(class_species_KO,', ')])

%%
% Unfortunately, the low number of L3, FS/BT cells does not allow for
% further analyses. Therefor we will focus on L5 pyramidal cells

Data_WT_L5 = Data(:,[2:5,8:17]);
Data_WT_L5 = cellfun(@(x) x(ind_mGluR3WT_L5_pp'), Data_WT_L5, 'UniformOutput', false);
Data_KO_L5 = Data(:,[2:5,8:17]);
Data_KO_L5 = cellfun(@(x) x(ind_mGluR3KO_L5_nn'), Data_KO_L5, 'UniformOutput', false);

%%
% We generate the model parameter histograms to eyeball the sample
% distributions and to get a first impression of differences between
% WT and KO.

figure()
for i=1:10
    subplot(2,5,i)
    EMin = min([Data_WT_L5{1,4+i},Data_KO_L5{1,4+i}]);
    EMax = max([Data_WT_L5{1,4+i},Data_KO_L5{1,4+i}]);
    edges = linspace(EMin,EMax,8);
    histogram(Data_WT_L5{1,4+i},edges,'FaceColor','b','Normalization','pdf')
    hold on;
    histogram(Data_KO_L5{1,4+i},edges,'FaceColor','r','Normalization','pdf')
    [f1,x1] = ksdensity(Data_WT_L5{1,4+i});
    [f2,x2] = ksdensity(Data_KO_L5{1,4+i});
    plot(x1,f1,'b','LineWidth',2)
    plot(x2,f2,'r','LineWidth',2)
    title(ColNam(4+i))
    xlim([EMin,EMax])
end
set(gcf,'units','centimeters','position',[1,1,25,15])
snapnow, close

%% 
% By eyeball inspection, most distributions seem to be nonnormal. Also, the
% distributions of WT and KO conditions very often look similar (even 
% though this is hard to tell given the low number of datat points). We cannot
% expect hugh differences in mean, median or standard deviation. This can
% also be seen in a boxplot summary:

figure()
for i=1:10
    subplot(2,5,i)
    Val = [Data_WT_L5{1,4+i},Data_KO_L5{1,4+i}];
    Label = [ones(1,length(Data_WT_L5{1,4+i})),2*ones(1,length(Data_KO_L5{1,4+i}))];
    boxplot(Val,Label,'Labels',{'WT','KO'},'Whisker',10)%,'Notch','on')hold on
    hold on;
    plot([mean(Val(Label==1)),mean(Val(Label==2))], 'xr')
    title(ColNam(4+i))
end
set(gcf,'units','centimeters','position',[1,1,25,15])
snapnow; close;

%% 
% We need to check if the mean and variance of the WT and KO
% distributions are significantly different. We do so by running a
% two-sample t-test (in case the distributions are halfway normal) and a
% Wilcoxon rank sum test (for non-normal distributions) for the mean, and
% a Bartlett’s, Brown-Forsythe and Levene’s test for the variance.

H_mean_t = zeros(1,10);
H_mean_rank = zeros(1,10);
H_var_brown = zeros(1,10);
H_var_levene = zeros(1,10);
H_var_bart = zeros(1,10);
for i=1:10
    x = Data_WT_L5{1,4+i};
    y = Data_KO_L5{1,4+i};
    % Check mean ...
    [h,~,ci,~] = ttest2(x,y,'VarType','unequal'); 
    H_mean_t(i) = h;
    [~,h] = ranksum(x,y);
    H_mean_rank(i) = h;
    % Check variance ...
    [p,~] = vartestn([x,y]',[zeros(1,length(x)),ones(1,length(y))]','TestType','BrownForsythe','Display','off');
    if p<0.05
        H_var_brown(i) = 1;
    end
    [p,~] = vartestn([x,y]',[zeros(1,length(x)),ones(1,length(y))]','TestType','LeveneQuadratic','Display','off');
    if p<0.05
        H_var_levene(i) = 1;
    end
    [p,~] = vartestn([x,y]',[zeros(1,length(x)),ones(1,length(y))]','TestType','Bartlett','Display','off');
    if p<0.05
        H_var_bart(i) = 1;
    end
end

id1 = find(H_mean_t==1);
id2 = find(H_mean_rank==1);

if (isempty(id1) && isempty(id2))
    disp('No mean is significantly different')
else
    for i=1:length(id1)
        disp(join(['According to a two-sample t-test, the mean of the' ColNam(4+id1(i)) 'distribution is significantly different between WT and KO']))
    end
    for i=1:length(id2)
        disp(join(['According to a Wilcoxon rank sum test, the mean of the' ColNam(4+id2(i)) 'distribution is significantly different between WT and KO']))
    end
end

id1 = find(H_var_bart==1);
id2 = find(H_var_brown==1);
id3 = find(H_var_levene==1);

if (isempty(id1) && isempty(id2))
    disp('No varaince is significantly different')
else
    for i=1:length(id1)
        disp(join(['According to a Bartlett’s test, the variance of the' ColNam(4+id1(i)) 'distribution is significantly different between WT and KO']))
    end
    for i=1:length(id2)
        disp(join(['According to a Brown-Forsythe test, the variance of the' ColNam(4+id2(i)) 'distribution is significantly different between WT and KO']))
    end
    for i=1:length(id3)
        disp(join(['According to a Levene’s test, the variance of the' ColNam(4+id3(i)) 'distribution is significantly different between WT and KO']))
    end
end

%%
% As some of the parameters are not independent, and in general the
% distribution can be regarded as multivariate, we perform a MANOVA test
% (under the probably wrong assumptions that the data comes from a n-dim
% normal distribution). As we have two groups (WT and KO),
% MANOVA is equivalent and the test reduces to Hotelling's T-square.

x = []; y = [];
for i=[1,2,4,5,6,7,9,10]
    x = [x Data_WT_L5{1,4+i}'];
    y = [y Data_KO_L5{1,4+i}'];   
end

z = [x;y];
group = [zeros(1,size(x,1)),ones(1,size(y,1))];
[~,p,~] = manova1(z,group);
[~,P] = HotellingsT2(x,y,0.05,0);

if (p<0.05 && P<0.05)
    disp("According to the MANOVA/Hotellings T2 test, the vector sample means are significantly different.")
elseif (p>=0.05 && P<0.05)
    disp("According to the MANOVA test, the vector of means of WT and KO conditions are not significantly different.")
elseif (p<0.05 && P>=0.05)
    disp("According to the Hotellings T2 test, the vector of means of WT and KO conditions are not significantly different.")
else
    disp("According to the MANOVA/Hotellings T2 test, the vector of means of WT and KO conditions are not significantly different.")
end


%% Analysis of WT (control) vs. KO (control) data
% In a first step we look at the number of data points with respect to
% specific conditions (age, cell type, layer, species)

N_WT = size(Data_WT_cont{1,1},1);
disp(['The WT data set contains ' num2str(N_WT) ' data points'])
class_age_WT = unique(Data_WT_cont{1,1});
disp(['The following age categories are present for WT: ' strjoin(class_age_WT,', ')])
class_celltype_WT = unique(Data_WT_cont{1,2});
disp(['The following cell types are present for control: ' strjoin(class_celltype_WT,', ')])
disp(['BB: ' num2str(length(ind_mGluR3WT_BB_Glu_control)) ', L3-PC: ' num2str(length(ind_mGluR3WT_L3_PC_Glu_control)) ', L5-PC: ' num2str(length(ind_mGluR3WT_L5_pp_Glu_control))])
class_layer_WT = unique(Data_WT_cont{1,3});
disp(['The following layers are present for WT: ', strjoin(class_layer_WT,', ')])
class_species_WT = unique(Data_WT_cont{1,4});
disp(['The recordings come from the follwoing species: ', strjoin(class_species_WT,', ')])
disp([' '])
N_KO = size(Data_KO_cont{1,1},1);
disp(['The KO data set contains ' num2str(N_KO) ' data points'])
class_age_KO = unique(Data_KO_cont{1,1});
disp(['The following age categories are present for KO: ' strjoin(class_age_KO,', ')])
class_celltype_KO = unique(Data_KO_cont{1,2});
disp(['The following cell types are present for KO: ' strjoin(class_celltype_KO,', ')])
disp(['FS: ' num2str(length(ind_mGluR3KO_FS_Glu_control)) ', L3-PC: ' num2str(length(ind_mGluR3KO_L3_PC_Glu_control)) ', L5-PC: ' num2str(length(ind_mGluR3KO_L5_nn_Glu_control))])
class_layer_KO = unique(Data_KO_cont{1,3});
disp(['The following layers are present for KO: ', strjoin(class_layer_KO,', ')])
class_species_KO = unique(Data_KO_cont{1,4});
disp(['The recordings come from the follwoing species: ', strjoin(class_species_KO,', ')])

%%
% Unfortunately, the low number of L3, FS/BT cells does not allow for
% further analyses. Therefor we will focus on L5 pyramidal cells

Data_WT_L5 = Data(:,[2:5,8:17]);
Data_WT_L5 = cellfun(@(x) x(ind_mGluR3WT_L5_pp_Glu_control'), Data_WT_L5, 'UniformOutput', false);
Data_KO_L5 = Data(:,[2:5,8:17]);
Data_KO_L5 = cellfun(@(x) x(ind_mGluR3KO_L5_nn_Glu_control'), Data_KO_L5, 'UniformOutput', false);

%%
% We generate the model parameter histograms to eyeball the sample
% distributions and to get a first impression of differences between
% WT and KO, and summarize the mean, median and std in a boxplot:

figure()
for i=1:10
    subplot(2,5,i)
    EMin = min([Data_WT_L5{1,4+i},Data_KO_L5{1,4+i}]);
    EMax = max([Data_WT_L5{1,4+i},Data_KO_L5{1,4+i}]);
    edges = linspace(EMin,EMax,8);
    histogram(Data_WT_L5{1,4+i},edges,'FaceColor','b','Normalization','pdf')
    hold on;
    histogram(Data_KO_L5{1,4+i},edges,'FaceColor','r','Normalization','pdf')
    [f1,x1] = ksdensity(Data_WT_L5{1,4+i});
    [f2,x2] = ksdensity(Data_KO_L5{1,4+i});
    plot(x1,f1,'b','LineWidth',2)
    plot(x2,f2,'r','LineWidth',2)
    title(ColNam(4+i))
    xlim([EMin,EMax])
end
set(gcf,'units','centimeters','position',[1,1,25,15])
snapnow, close

figure()
for i=1:10
    subplot(2,5,i)
    Val = [Data_WT_L5{1,4+i},Data_KO_L5{1,4+i}];
    Label = [ones(1,length(Data_WT_L5{1,4+i})),2*ones(1,length(Data_KO_L5{1,4+i}))];
    boxplot(Val,Label,'Labels',{'WT','KO'},'Whisker',10)%,'Notch','on')hold on
    hold on;
    plot([mean(Val(Label==1)),mean(Val(Label==2))], 'xr')
    title(ColNam(4+i))
end
set(gcf,'units','centimeters','position',[1,1,25,15])
snapnow; close;

%% 
% We need to check if the mean and variance of the WT and KO
% distributions are significantly different. We do so by running a
% two-sample t-test (in case the distributions are halfway normal) and a
% Wilcoxon rank sum test (for non-normal distributions) for the mean, and
% a Bartlett’s, Brown-Forsythe and Levene’s test for the variance.

H_mean_t = zeros(1,10);
H_mean_rank = zeros(1,10);
H_var_brown = zeros(1,10);
H_var_levene = zeros(1,10);
H_var_bart = zeros(1,10);
for i=1:10
    x = Data_WT_L5{1,4+i};
    y = Data_KO_L5{1,4+i};
    % Check mean ...
    [h,~,ci,~] = ttest2(x,y,'VarType','unequal'); 
    H_mean_t(i) = h;
    [~,h] = ranksum(x,y);
    H_mean_rank(i) = h;
    % Check variance ...
    [p,~] = vartestn([x,y]',[zeros(1,length(x)),ones(1,length(y))]','TestType','BrownForsythe','Display','off');
    if p<0.05
        H_var_brown(i) = 1;
    end
    [p,~] = vartestn([x,y]',[zeros(1,length(x)),ones(1,length(y))]','TestType','LeveneQuadratic','Display','off');
    if p<0.05
        H_var_levene(i) = 1;
    end
    [p,~] = vartestn([x,y]',[zeros(1,length(x)),ones(1,length(y))]','TestType','Bartlett','Display','off');
    if p<0.05
        H_var_bart(i) = 1;
    end
end

id1 = find(H_mean_t==1);
id2 = find(H_mean_rank==1);

if (isempty(id1) && isempty(id2))
    disp('No mean is significantly different')
else
    for i=1:length(id1)
        disp(join(['According to a two-sample t-test, the mean of the' ColNam(4+id1(i)) 'distribution is significantly different between WT and KO']))
    end
    for i=1:length(id2)
        disp(join(['According to a Wilcoxon rank sum test, the mean of the' ColNam(4+id2(i)) 'distribution is significantly different between WT and KO']))
    end
end

id1 = find(H_var_bart==1);
id2 = find(H_var_brown==1);
id3 = find(H_var_levene==1);

if (isempty(id1) && isempty(id2))
    disp('No varaince is significantly different')
else
    for i=1:length(id1)
        disp(join(['According to a Bartlett’s test, the variance of the' ColNam(4+id1(i)) 'distribution is significantly different between WT and KO']))
    end
    for i=1:length(id2)
        disp(join(['According to a Brown-Forsythe test, the variance of the' ColNam(4+id2(i)) 'distribution is significantly different between WT and KO']))
    end
    for i=1:length(id3)
        disp(join(['According to a Levene’s test, the variance of the' ColNam(4+id3(i)) 'distribution is significantly different between WT and KO']))
    end
end

%%
% As some of the parameters are not independent, and in general the
% distribution can be regarded as multivariate, we perform a MANOVA test
% (under the probably wrong assumptions that the data comes from a n-dim
% normal distribution). As we have two groups (WT and KO),
% MANOVA is equivalent and the test reduces to Hotelling's T-square.

x = []; y = [];
for i=[1,2,4,5,6,7,9,10]
    x = [x Data_WT_L5{1,4+i}'];
    y = [y Data_KO_L5{1,4+i}'];   
end

z = [x;y];
group = [zeros(1,size(x,1)),ones(1,size(y,1))];
[~,p,~] = manova1(z,group);
[~,P] = HotellingsT2(x,y,0.05,0);

if (p<0.05 && P<0.05)
    disp("According to the MANOVA/Hotellings T2 test, the vector sample means are significantly different.")
elseif (p>=0.05 && P<0.05)
    disp("According to the MANOVA test, the vector of means of WT and KO conditions are not significantly different.")
elseif (p<0.05 && P>=0.05)
    disp("According to the Hotellings T2 test, the vector of means of WT and KO conditions are not significantly different.")
else
    disp("According to the MANOVA/Hotellings T2 test, the vector of means of WT and KO conditions are not significantly different.")
end

%% Analysis of WT (agonist) vs. KO (agonist) data
% In a first step we look at the number of data points with respect to
% specific conditions (age, cell type, layer, species)

N_WT = size(Data_WT_agon{1,1},1);
disp(['The WT data set contains ' num2str(N_WT) ' data points'])
class_age_WT = unique(Data_WT_agon{1,1});
disp(['The following age categories are present for WT: ' strjoin(class_age_WT,', ')])
class_celltype_WT = unique(Data_WT_agon{1,2});
disp(['The following cell types are present for control: ' strjoin(class_celltype_WT,', ')])
disp(['BB: ' num2str(length(ind_mGluR3WT_BB_Glu_agonist)) ', L3-PC: ' num2str(length(ind_mGluR3WT_L3_PC_Glu_agonist)) ', L5-PC: ' num2str(length(ind_mGluR3WT_L5_pp_Glu_agonist))])
class_layer_WT = unique(Data_WT_agon{1,3});
disp(['The following layers are present for WT: ', strjoin(class_layer_WT,', ')])
class_species_WT = unique(Data_WT_agon{1,4});
disp(['The recordings come from the follwoing species: ', strjoin(class_species_WT,', ')])
disp([' '])
N_KO = size(Data_KO_agon{1,1},1);
disp(['The KO data set contains ' num2str(N_KO) ' data points'])
class_age_KO = unique(Data_KO_agon{1,1});
disp(['The following age categories are present for KO: ' strjoin(class_age_KO,', ')])
class_celltype_KO = unique(Data_KO_agon{1,2});
disp(['The following cell types are present for KO: ' strjoin(class_celltype_KO,', ')])
disp(['FS: ' num2str(length(ind_mGluR3KO_FS_Glu_agonist)) ', L3-PC: ' num2str(length(ind_mGluR3KO_L3_PC_Glu_agonist)) ', L5-PC: ' num2str(length(ind_mGluR3KO_L5_nn_Glu_agonist))])
class_layer_KO = unique(Data_KO_agon{1,3});
disp(['The following layers are present for KO: ', strjoin(class_layer_KO,', ')])
class_species_KO = unique(Data_KO_agon{1,4});
disp(['The recordings come from the follwoing species: ', strjoin(class_species_KO,', ')])

%%
% Unfortunately, the low number of L3, FS/BT cells does not allow for
% further analyses. Therefor we will focus on L5 pyramidal cells

Data_WT_L5 = Data(:,[2:5,8:17]);
Data_WT_L5 = cellfun(@(x) x(ind_mGluR3WT_L5_pp_Glu_agonist'), Data_WT_L5, 'UniformOutput', false);
Data_KO_L5 = Data(:,[2:5,8:17]);
Data_KO_L5 = cellfun(@(x) x(ind_mGluR3KO_L5_nn_Glu_agonist'), Data_KO_L5, 'UniformOutput', false);

%%
% We generate the model parameter histograms to eyeball the sample
% distributions and to get a first impression of differences between
% WT and KO, and summarize the mean, median and std in a boxplot:

figure()
for i=1:10
    subplot(2,5,i)
    EMin = min([Data_WT_L5{1,4+i},Data_KO_L5{1,4+i}]);
    EMax = max([Data_WT_L5{1,4+i},Data_KO_L5{1,4+i}]);
    edges = linspace(EMin,EMax,8);
    histogram(Data_WT_L5{1,4+i},edges,'FaceColor','b','Normalization','pdf')
    hold on;
    histogram(Data_KO_L5{1,4+i},edges,'FaceColor','r','Normalization','pdf')
    [f1,x1] = ksdensity(Data_WT_L5{1,4+i});
    [f2,x2] = ksdensity(Data_KO_L5{1,4+i});
    plot(x1,f1,'b','LineWidth',2)
    plot(x2,f2,'r','LineWidth',2)
    title(ColNam(4+i))
    xlim([EMin,EMax])
end
set(gcf,'units','centimeters','position',[1,1,25,15])
snapnow, close

figure()
for i=1:10
    subplot(2,5,i)
    Val = [Data_WT_L5{1,4+i},Data_KO_L5{1,4+i}];
    Label = [ones(1,length(Data_WT_L5{1,4+i})),2*ones(1,length(Data_KO_L5{1,4+i}))];
    boxplot(Val,Label,'Labels',{'WT','KO'},'Whisker',10)%,'Notch','on')hold on
    hold on;
    plot([mean(Val(Label==1)),mean(Val(Label==2))], 'xr')
    title(ColNam(4+i))
end
set(gcf,'units','centimeters','position',[1,1,25,15])
snapnow; close;

%% 
% We need to check if the mean and variance of the WT and KO
% distributions are significantly different. We do so by running a
% two-sample t-test (in case the distributions are halfway normal) and a
% Wilcoxon rank sum test (for non-normal distributions) for the mean, and
% a Bartlett’s, Brown-Forsythe and Levene’s test for the variance.

H_mean_t = zeros(1,10);
H_mean_rank = zeros(1,10);
H_var_brown = zeros(1,10);
H_var_levene = zeros(1,10);
H_var_bart = zeros(1,10);
for i=1:10
    x = Data_WT_L5{1,4+i};
    y = Data_KO_L5{1,4+i};
    % Check mean ...
    [h,~,ci,~] = ttest2(x,y,'VarType','unequal'); 
    H_mean_t(i) = h;
    [~,h] = ranksum(x,y);
    H_mean_rank(i) = h;
    % Check variance ...
    [p,~] = vartestn([x,y]',[zeros(1,length(x)),ones(1,length(y))]','TestType','BrownForsythe','Display','off');
    if p<0.05
        H_var_brown(i) = 1;
    end
    [p,~] = vartestn([x,y]',[zeros(1,length(x)),ones(1,length(y))]','TestType','LeveneQuadratic','Display','off');
    if p<0.05
        H_var_levene(i) = 1;
    end
    [p,~] = vartestn([x,y]',[zeros(1,length(x)),ones(1,length(y))]','TestType','Bartlett','Display','off');
    if p<0.05
        H_var_bart(i) = 1;
    end
end

id1 = find(H_mean_t==1);
id2 = find(H_mean_rank==1);

if (isempty(id1) && isempty(id2))
    disp('No mean is significantly different')
else
    for i=1:length(id1)
        disp(join(['According to a two-sample t-test, the mean of the' ColNam(4+id1(i)) 'distribution is significantly different between WT and KO']))
    end
    for i=1:length(id2)
        disp(join(['According to a Wilcoxon rank sum test, the mean of the' ColNam(4+id2(i)) 'distribution is significantly different between WT and KO']))
    end
end

id1 = find(H_var_bart==1);
id2 = find(H_var_brown==1);
id3 = find(H_var_levene==1);

if (isempty(id1) && isempty(id2))
    disp('No varaince is significantly different')
else
    for i=1:length(id1)
        disp(join(['According to a Bartlett’s test, the variance of the' ColNam(4+id1(i)) 'distribution is significantly different between WT and KO']))
    end
    for i=1:length(id2)
        disp(join(['According to a Brown-Forsythe test, the variance of the' ColNam(4+id2(i)) 'distribution is significantly different between WT and KO']))
    end
    for i=1:length(id3)
        disp(join(['According to a Levene’s test, the variance of the' ColNam(4+id3(i)) 'distribution is significantly different between WT and KO']))
    end
end

%%
% As some of the parameters are not independent, and in general the
% distribution can be regarded as multivariate, we perform a MANOVA test
% (under the probably wrong assumptions that the data comes from a n-dim
% normal distribution). As we have two groups (WT and KO),
% MANOVA is equivalent and the test reduces to Hotelling's T-square.

x = []; y = [];
for i=[1,2,4,5,6,7,9,10]
    x = [x Data_WT_L5{1,4+i}'];
    y = [y Data_KO_L5{1,4+i}'];   
end

z = [x;y];
group = [zeros(1,size(x,1)),ones(1,size(y,1))];
[~,p,~] = manova1(z,group);
[~,P] = HotellingsT2(x,y,0.05,0);

if (p<0.05 && P<0.05)
    disp("According to the MANOVA/Hotellings T2 test, the vector sample means are significantly different.")
elseif (p>=0.05 && P<0.05)
    disp("According to the MANOVA test, the vector of means of WT and KO conditions are not significantly different.")
elseif (p<0.05 && P>=0.05)
    disp("According to the Hotellings T2 test, the vector of means of WT and KO conditions are not significantly different.")
else
    disp("According to the MANOVA/Hotellings T2 test, the vector of means of WT and KO conditions are not significantly different.")
end

%% Analysis of WT (control) vs. WT (agonist) data
% In a first step we look at the number of data points with respect to
% specific conditions (age, cell type, layer, species)

N_WT_cont = size(Data_WT_cont{1,1},1);
disp(['The WT/control data set contains ' num2str(N_WT_cont) ' data points'])
class_age_WT_cont = unique(Data_WT_cont{1,1});
disp(['The following age categories are present for WT/control: ' strjoin(class_age_WT_cont,', ')])
class_celltype_WT_cont = unique(Data_WT_cont{1,2});
disp(['The following cell types are present for WT/control: ' strjoin(class_celltype_WT_cont,', ')])
disp(['BB: ' num2str(length(ind_mGluR3WT_BB_Glu_control)) ', L3-PC: ' num2str(length(ind_mGluR3WT_L3_PC_Glu_control)) ', L5-PC: ' num2str(length(ind_mGluR3WT_L5_pp_Glu_control))])
class_layer_WT_cont = unique(Data_WT_cont{1,3});
disp(['The following layers are present for WT/control: ', strjoin(class_layer_WT_cont,', ')])
class_species_WT_cont = unique(Data_WT_cont{1,4});
disp(['The recordings come from the follwoing species: ', strjoin(class_species_WT_cont,', ')])
disp([' '])
N_WT_agon = size(Data_WT_agon{1,1},1);
disp(['The WT/agonist data set contains ' num2str(N_WT_agon) ' data points'])
class_age_WT_agon = unique(Data_WT_agon{1,1});
disp(['The following age categories are present for WT/agonist: ' strjoin(class_age_WT_agon,', ')])
class_celltype_WT_agon = unique(Data_WT_agon{1,2});
disp(['The following cell types are present for WT/agon: ' strjoin(class_celltype_WT_agon,', ')])
disp(['BB: ' num2str(length(ind_mGluR3WT_BB_Glu_agonist)) ', L3-PC: ' num2str(length(ind_mGluR3WT_L3_PC_Glu_agonist)) ', L5-PC: ' num2str(length(ind_mGluR3WT_L5_pp_Glu_agonist))])
class_layer_WT_agon = unique(Data_WT_agon{1,3});
disp(['The following layers are present for WT/agonist: ', strjoin(class_layer_WT_agon,', ')])
class_species_WT_agon = unique(Data_WT_agon{1,4});
disp(['The recordings come from the follwoing species: ', strjoin(class_species_WT_agon,', ')])

%% Analysis of KO (control) vs. KO (agonist) data
% In a first step we look at the number of data points with respect to
% specific conditions (age, cell type, layer, species)

N_KO_cont = size(Data_KO_cont{1,1},1);
disp(['The KO/control data set contains ' num2str(N_KO_cont) ' data points'])
class_age_KO_cont = unique(Data_KO_cont{1,1});
disp(['The following age categories are present for KO/control: ' strjoin(class_age_KO_cont,', ')])
class_celltype_KO_cont = unique(Data_KO_cont{1,2});
disp(['The following cell types are present for KO/control: ' strjoin(class_celltype_KO_cont,', ')])
disp(['FS: ' num2str(length(ind_mGluR3KO_FS_Glu_control)) ', L3-PC: ' num2str(length(ind_mGluR3KO_L3_PC_Glu_control)) ', L5-PC: ' num2str(length(ind_mGluR3KO_L5_nn_Glu_control))])
class_layer_KO_cont = unique(Data_KO_cont{1,3});
disp(['The following layers are present for KO/control: ', strjoin(class_layer_KO_cont,', ')])
class_species_KO_cont = unique(Data_KO_cont{1,4});
disp(['The recordings come from the follwoing species: ', strjoin(class_species_KO_cont,', ')])
disp([' '])
N_KO_agon = size(Data_KO_agon{1,1},1);
disp(['The KO/agonist data set contains ' num2str(N_KO_agon) ' data points'])
class_age_KO_agon = unique(Data_KO_agon{1,1});
disp(['The following age categories are present for KO/agonist: ' strjoin(class_age_KO_agon,', ')])
class_celltype_KO_agon = unique(Data_KO_agon{1,2});
disp(['The following cell types are present for KO/agon: ' strjoin(class_celltype_KO_agon,', ')])
disp(['FS: ' num2str(length(ind_mGluR3KO_FS_Glu_agonist)) ', L3-PC: ' num2str(length(ind_mGluR3KO_L3_PC_Glu_agonist)) ', L5-PC: ' num2str(length(ind_mGluR3KO_L5_nn_Glu_agonist))])
class_layer_KO_agon = unique(Data_KO_agon{1,3});
disp(['The following layers are present for KO/agonist: ', strjoin(class_layer_KO_agon,', ')])
class_species_KO_agon = unique(Data_KO_agon{1,4});
disp(['The recordings come from the follwoing species: ', strjoin(class_species_KO_agon,', ')])
