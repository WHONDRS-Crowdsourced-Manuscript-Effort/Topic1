% visualize some of the aggregate properties for core and satellite species
% CM, Mar 11, 2022

%% Import data from text file

% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 51);
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["VarName1", "id", "Mass", "MolForm", "C", "H", "O", "N", "C13", "S", "P", "Na", "El_comp", "Class", "NeutralMass", "Error_ppm", "Candidates", "AI", "AI_Mod", "DBE", "DBE_O", "DBE_AI", "GFE", "kmassCH2", "kdefectCH2", "NOSC", "OtoC_ratio", "HtoC_ratio", "NtoC_ratio", "PtoC_ratio", "NtoP_ratio", "bs1_class", "bs2_class", "bs3_class", "delGcox0PerCmol", "delGcoxPerCmol", "lamO20", "lamO2", "delGd0", "delGd", "nmf", "occupancy_sed", "occupancy_water", "percoccup_sed", "percoccup_water", "csflagemergent_sed", "csflagemergent_water", "csflagpca_sed", "csflagpca_water", "csflagrf_sed", "csflagrf_water"];
opts.VariableTypes = ["double", "double", "double", "string", "double", "double", "double", "double", "double", "double", "double", "double", "categorical", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "categorical", "categorical", "categorical", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts = setvaropts(opts, "MolForm", "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["MolForm", "El_comp", "Class", "NtoP_ratio", "bs1_class", "bs2_class", "bs3_class", "csflagemergent_sed", "csflagemergent_water", "csflagpca_sed", "csflagpca_water", "csflagrf_sed", "csflagrf_water"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, "id", "TrimNonNumeric", true);
opts = setvaropts(opts, "id", "ThousandsSeparator", ",");

% select rep1 or rep2 - this is data downloaded from github, https://github.com/WHONDRS-Crowdsourced-Manuscript-Effort/Topic1/tree/main/4_gather.thresholds
data = readtable("FTICR_crosstable_rep.merged1_all_em.thres_2022-03-07.csv", opts);
% data = readtable("FTICR_crosstable_rep.merged2_all_em.thres_2022-03-07.csv", opts);


%% visualize the data with violin or swarm plots

% emergent
cat1 = ["In-between" "Satellite" "Core"];
%  water
x1 = categorical(data.csflagemergent_water,cat1);
figure1 = figure;
subplot(7,1,1), swarmchart(x1,data.DBE,'.'), title('DBE')
subplot(7,1,2), swarmchart(x1,data.DBE_O,'.'), title('DBE O')
subplot(7,1,3), swarmchart(x1,data.DBE_AI,'.'), title('DBE AI')
subplot(7,1,4), swarmchart(x1,data.AI,'.'), title('AI')
subplot(7,1,5), swarmchart(x1,data.AI_Mod,'.'), title('AI mod')
subplot(7,1,6), swarmchart(x1,data.GFE,'.'), title('GFE')
subplot(7,1,7), swarmchart(x1,data.NOSC,'.'), title('NOSC'), hold on
annotation(figure1,'textbox', [0.13 0.94 0.16 0.02],'String',{'Water, emergent'},'FitBoxToText','on');
%  sediment
x1 = categorical(data.csflagemergent_sed,cat1);
figure1 = figure;
subplot(7,1,1), swarmchart(x1,data.DBE,'.'), title('DBE')
subplot(7,1,2), swarmchart(x1,data.DBE_O,'.'), title('DBE O')
subplot(7,1,3), swarmchart(x1,data.DBE_AI,'.'), title('DBE AI')
subplot(7,1,4), swarmchart(x1,data.AI,'.'), title('AI')
subplot(7,1,5), swarmchart(x1,data.AI_Mod,'.'), title('AI mod')
subplot(7,1,6), swarmchart(x1,data.GFE,'.'), title('GFE')
subplot(7,1,7), swarmchart(x1,data.NOSC,'.'), title('NOSC'), hold on
annotation(figure1,'textbox', [0.13 0.94 0.16 0.02],'String',{'Sediment, emergent'},'FitBoxToText','on');

% pca
cat1 = ["Satellite" "Core"];
%  water
x1 = categorical(data.csflagpca_water,cat1);
figure1 = figure;
subplot(7,1,1), swarmchart(x1,data.DBE,'.'), title('DBE')
subplot(7,1,2), swarmchart(x1,data.DBE_O,'.'), title('DBE O')
subplot(7,1,3), swarmchart(x1,data.DBE_AI,'.'), title('DBE AI')
subplot(7,1,4), swarmchart(x1,data.AI,'.'), title('AI')
subplot(7,1,5), swarmchart(x1,data.AI_Mod,'.'), title('AI mod')
subplot(7,1,6), swarmchart(x1,data.GFE,'.'), title('GFE')
subplot(7,1,7), swarmchart(x1,data.NOSC,'.'), title('NOSC'), hold on
annotation(figure1,'textbox', [0.13 0.94 0.16 0.02],'String',{'Water, PCA'},'FitBoxToText','on');
%  sediment
x1 = categorical(data.csflagpca_sed,cat1);
figure1 = figure;
subplot(7,1,1), swarmchart(x1,data.DBE,'.'), title('DBE')
subplot(7,1,2), swarmchart(x1,data.DBE_O,'.'), title('DBE O')
subplot(7,1,3), swarmchart(x1,data.DBE_AI,'.'), title('DBE AI')
subplot(7,1,4), swarmchart(x1,data.AI,'.'), title('AI')
subplot(7,1,5), swarmchart(x1,data.AI_Mod,'.'), title('AI mod')
subplot(7,1,6), swarmchart(x1,data.GFE,'.'), title('GFE')
subplot(7,1,7), swarmchart(x1,data.NOSC,'.'), title('NOSC'), hold on
annotation(figure1,'textbox', [0.13 0.94 0.16 0.02],'String',{'Sediment, PCA'},'FitBoxToText','on');

% randomforest
cat1 = ["Satellite" "Core"];
%  water
x1 = categorical(data.csflagrf_water,cat1);
figure1 = figure;
subplot(7,1,1), swarmchart(x1,data.DBE,'.'), title('DBE')
subplot(7,1,2), swarmchart(x1,data.DBE_O,'.'), title('DBE O')
subplot(7,1,3), swarmchart(x1,data.DBE_AI,'.'), title('DBE AI')
subplot(7,1,4), swarmchart(x1,data.AI,'.'), title('AI')
subplot(7,1,5), swarmchart(x1,data.AI_Mod,'.'), title('AI mod')
subplot(7,1,6), swarmchart(x1,data.GFE,'.'), title('GFE')
subplot(7,1,7), swarmchart(x1,data.NOSC,'.'), title('NOSC'), hold on
annotation(figure1,'textbox', [0.13 0.94 0.16 0.02],'String',{'Water, Random Forest'},'FitBoxToText','on');
%  sediment
x1 = categorical(data.csflagrf_sed,cat1);
figure1 = figure;
subplot(7,1,1), swarmchart(x1,data.DBE,'.'), title('DBE')
subplot(7,1,2), swarmchart(x1,data.DBE_O,'.'), title('DBE O')
subplot(7,1,3), swarmchart(x1,data.DBE_AI,'.'), title('DBE AI')
subplot(7,1,4), swarmchart(x1,data.AI,'.'), title('AI')
subplot(7,1,5), swarmchart(x1,data.AI_Mod,'.'), title('AI mod')
subplot(7,1,6), swarmchart(x1,data.GFE,'.'), title('GFE')
subplot(7,1,7), swarmchart(x1,data.NOSC,'.'), title('NOSC'), hold on
annotation(figure1,'textbox', [0.13 0.94 0.16 0.02],'String',{'Sediment, Random Forest'},'FitBoxToText','on');


%% NMDS: too large
% X = [data.AI; data.AI_Mod; data.DBE; data.DBE_O; data.DBE_AI; data.GFE; data.NOSC; data.kdefectCH2];
% dissimilarities = pdist(X);
% [Y,stress,disparities] = mdscale(dissimilarities,2);
% distances = pdist(Y);
% [dum,ord] = sortrows([disparities(:) dissimilarities(:)]);
% plot(dissimilarities,distances,'bo', ...
% dissimilarities(ord),disparities(ord),'r.-');
% xlabel('Dissimilarities'); ylabel('Distances/Disparities')
% legend({'Distances' 'Disparities'},'Location','NW');

%% cluster: not very useful
% yy=[data.DBE_O,data.DBE_AI,data.AI,data.AI_Mod,data.GFE,data.NOSC];
% idx = kmeans(yy,2);
% dist = pdist(yy);
% Z = linkage(dist);
% dendrogram(Z)

%% PCA
% run a PCA on ALL the data together, then plot the different subgroups of the data
% - base grouping: sediment core, sediment satellite, water core, water satellite
% - additional: plot MF that are both in core water & core sediment, core water & satellite sediment, 
%                    core sediment & satellite water, satellite water & satellite sediment
X = [data.AI, data.AI_Mod, data.DBE, data.DBE_O, data.DBE_AI, data.GFE, data.NOSC, data.kdefectCH2];
vbls={'AI','AI_Mod','DBE','DBE_O','DBE_AI','GFE','NOSC','kdefectCH2'}
Z = zscore(X); % Standardized data
[coefs,score] = pca(Z);
ndim = 2; % dimensions of the biplot
% emergent
cat1 = ["In-between" "Satellite" "Core"];
figure, 
x1 = categorical(data.csflagemergent_water,cat1);
subplot(3,2,1), biplot(coefs(:,1:ndim),'Scores',score(x1==cat1(1),1:ndim),'VarLabels',vbls); title('emergent,water,in-between')
subplot(3,2,3), biplot(coefs(:,1:ndim),'Scores',score(x1==cat1(2),1:ndim),'VarLabels',vbls); title('emergent,water,satellite')
subplot(3,2,5), biplot(coefs(:,1:ndim),'Scores',score(x1==cat1(3),1:ndim),'VarLabels',vbls); title('emergent,water,core')
x1 = categorical(data.csflagemergent_sed,cat1);
subplot(3,2,2), biplot(coefs(:,1:ndim),'Scores',score(x1==cat1(1),1:ndim),'VarLabels',vbls); title('emergent,sed,in-between')
subplot(3,2,4), biplot(coefs(:,1:ndim),'Scores',score(x1==cat1(2),1:ndim),'VarLabels',vbls); title('emergent,sed,satellite')
subplot(3,2,6), biplot(coefs(:,1:ndim),'Scores',score(x1==cat1(3),1:ndim),'VarLabels',vbls); title('emergent,sed,core')

x1 = categorical(data.csflagemergent_water,cat1);
x2 = categorical(data.csflagemergent_sed,cat1);
ind_cwcd = x1==cat1(3)&x2==cat1(3);
ind_cwsd = x1==cat1(3)&x2==cat1(2);
ind_swcd = x1==cat1(2)&x2==cat1(3);
ind_swsd = x1==cat1(2)&x2==cat1(2); 
figure, 
subplot(2,2,1), biplot(coefs(:,1:ndim),'Scores',score(ind_cwcd,1:ndim),'VarLabels',vbls); title('emergent, core water & core sed')
subplot(2,2,2), biplot(coefs(:,1:ndim),'Scores',score(ind_cwsd,1:ndim),'VarLabels',vbls); title('emergent, core water & sat sed')
subplot(2,2,3), biplot(coefs(:,1:ndim),'Scores',score(ind_swcd,1:ndim),'VarLabels',vbls); title('emergent, sat water & core sed')
subplot(2,2,4), biplot(coefs(:,1:ndim),'Scores',score(ind_swsd,1:ndim),'VarLabels',vbls); title('emergent, sat water & sat sed')


% pca
figure, 
cat1 = ["Satellite" "Core"];
x1 = categorical(data.csflagpca_water,cat1);
subplot(2,2,1), biplot(coefs(:,1:ndim),'Scores',score(x1==cat1(1),1:ndim),'VarLabels',vbls); title('pca,water,satellite')
subplot(2,2,3), biplot(coefs(:,1:ndim),'Scores',score(x1==cat1(2),1:ndim),'VarLabels',vbls);title('pca,water,core')
x1 = categorical(data.csflagemergent_sed,cat1);
subplot(2,2,2), biplot(coefs(:,1:ndim),'Scores',score(x1==cat1(1),1:ndim),'VarLabels',vbls); title('pca,sed,satellite')
subplot(2,2,4), biplot(coefs(:,1:ndim),'Scores',score(x1==cat1(2),1:ndim),'VarLabels',vbls);title('pca,sed,core')

x1 = categorical(data.csflagpca_water,cat1);
x2 = categorical(data.csflagpca_sed,cat1);
ind_cwcd = x1==cat1(2)&x2==cat1(2);
ind_cwsd = x1==cat1(2)&x2==cat1(1);
ind_swcd = x1==cat1(1)&x2==cat1(2);
ind_swsd = x1==cat1(1)&x2==cat1(1); 
figure, 
subplot(2,2,1), biplot(coefs(:,1:ndim),'Scores',score(ind_cwcd,1:ndim),'VarLabels',vbls); title('pca, core water & core sed')
subplot(2,2,2), biplot(coefs(:,1:ndim),'Scores',score(ind_cwsd,1:ndim),'VarLabels',vbls); title('pca, core water & sat sed')
subplot(2,2,3), biplot(coefs(:,1:ndim),'Scores',score(ind_swcd,1:ndim),'VarLabels',vbls); title('pca, sat water & core sed')
subplot(2,2,4), biplot(coefs(:,1:ndim),'Scores',score(ind_swsd,1:ndim),'VarLabels',vbls); title('pca, sat water & sat sed')


% rf
figure, 
cat1 = ["Satellite" "Core"];
x1 = categorical(data.csflagrf_water,cat1);
subplot(2,2,1), biplot(coefs(:,1:ndim),'Scores',score(x1==cat1(1),1:ndim),'VarLabels',vbls); title('rf,water,satellite')
subplot(2,2,3), biplot(coefs(:,1:ndim),'Scores',score(x1==cat1(2),1:ndim),'VarLabels',vbls);title('rf,water,core')
x1 = categorical(data.csflagrf_sed,cat1);
subplot(2,2,2), biplot(coefs(:,1:ndim),'Scores',score(x1==cat1(1),1:ndim),'VarLabels',vbls); title('rf,sed,satellite')
subplot(2,2,4), biplot(coefs(:,1:ndim),'Scores',score(x1==cat1(2),1:ndim),'VarLabels',vbls);title('rf,sed,core')

x1 = categorical(data.csflagrf_water,cat1);
x2 = categorical(data.csflagrf_sed,cat1);
ind_cwcd = x1==cat1(2)&x2==cat1(2);
ind_cwsd = x1==cat1(2)&x2==cat1(1);
ind_swcd = x1==cat1(1)&x2==cat1(2);
ind_swsd = x1==cat1(1)&x2==cat1(1); 
figure, 
subplot(2,2,1), biplot(coefs(:,1:ndim),'Scores',score(ind_cwcd,1:ndim),'VarLabels',vbls); title('rf, core water & core sed')
subplot(2,2,2), biplot(coefs(:,1:ndim),'Scores',score(ind_cwsd,1:ndim),'VarLabels',vbls); title('rf, core water & sat sed')
subplot(2,2,3), biplot(coefs(:,1:ndim),'Scores',score(ind_swcd,1:ndim),'VarLabels',vbls); title('rf, sat water & core sed')
subplot(2,2,4), biplot(coefs(:,1:ndim),'Scores',score(ind_swsd,1:ndim),'VarLabels',vbls); title('rf, sat water & sat sed')

%% supervised learning to separate core from satellite - to do

