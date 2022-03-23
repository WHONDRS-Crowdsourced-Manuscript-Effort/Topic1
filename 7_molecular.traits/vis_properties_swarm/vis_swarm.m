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


%% visualize

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

