% run a PCA 
% CM, Mar 23, 2022

%% Import data from text file
% https://github.com/WHONDRS-Crowdsourced-Manuscript-Effort/Topic1/tree/main/4_gather.thresholds
opts = delimitedTextImportOptions("NumVariables", 56);
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["X", "id", "Mass", "MolForm", "C", "H", "O", "N", "C13", "S", "P", "Na", "El_comp", "Class", "NeutralMass", "Error_ppm", "Candidates", "AI", "AI_Mod", "DBE", "DBE_O", "DBE_AI", "GFE", "kmassCH2", "kdefectCH2", "NOSC", "OtoC_ratio", "HtoC_ratio", "NtoC_ratio", "PtoC_ratio", "NtoP_ratio", "bs1_class", "bs2_class", "bs3_class", "delGcox0PerCmol", "delGcoxPerCmol", "lamO20", "lamO2", "delGd0", "delGd", "nmf", "occupancy_sed", "occupancy_water", "percoccup_sed", "percoccup_water", "csflagemergent_sed", "csflagemergent_water", "csflagpca_sed", "csflagpca_water", "csflagrf_sed", "csflagrf_water", "habitatoverlap", "csflagemergent_overlap", "csflagemergent_generaloverlap", "csflagpca_generaloverlap", "csflagrf_generaloverlap"];
opts.VariableTypes = ["double", "double", "double", "string", "double", "double", "double", "double", "double", "double", "double", "double", "categorical", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "categorical", "categorical", "categorical", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts = setvaropts(opts, "MolForm", "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["MolForm", "El_comp", "Class", "NtoP_ratio", "bs1_class", "bs2_class", "bs3_class", "csflagemergent_sed", "csflagemergent_water", "csflagpca_sed", "csflagpca_water", "csflagrf_sed", "csflagrf_water", "habitatoverlap", "csflagemergent_overlap", "csflagemergent_generaloverlap", "csflagpca_generaloverlap", "csflagrf_generaloverlap"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, "id", "TrimNonNumeric", true);
opts = setvaropts(opts, "id", "ThousandsSeparator", ",");


% data = readtable("FTICR_crosstable_rep.merged2_all_em.thres_2022-03-23.csv", opts);
data = readtable("FTICR_crosstable_rep.merged1_all_em.thres_2022-03-23.csv", opts);

% select method used to decide on Core/Satellite designation
method = 0; % 0: emergent, 1: pca, 2: random forest

switch method
    case 0      % emergent
        ti1 = 'emergent';
        cat1 = ["Satellite" "Core" "In-between"];
        x1w = categorical(data.csflagemergent_water,cat1);  
        x1s = categorical(data.csflagemergent_sed,cat1); 

        % core water and sediment, or satellite in both water and sediment
        cat2 = ["Global core" "Global satellite","Shifter"];
        x1c = categorical(data.csflagemergent_generaloverlap,cat2);  
        ind_cwcd = x1w==cat1(2)&x1s==cat1(2);
        ind_swsd = x1w==cat1(1)&x1s==cat1(1);
    case 1      % pca
        ti1 = 'pca';
        cat1 = ["Satellite" "Core"];
        x1w = categorical(data.csflagpca_water,cat1);  
        x1s = categorical(data.csflagpca_sed,cat1);    

        cat2 = ["Global core" "Global satellite" "Shifter"];
        x1c = categorical(data.csflagpca_generaloverlap,cat2);  
        ind_cwcd = x1w==cat1(2)&x1s==cat1(2); % ind_cwcd = (x1c==cat2(1));
        ind_swsd = x1w==cat1(1)&x1s==cat1(1); % ind_swsd = (x1c==cat2(2));
    case 2      % random forest
        ti1 = 'random forest';
        cat1 = ["Satellite" "Core"];
        x1w = categorical(data.csflagrf_water,cat1);       
        x1s = categorical(data.csflagrf_sed,cat1);       

        cat2 = ["Global core" "Global satellite" "Shifter"];
        x1c = categorical(data.csflagrf_generaloverlap,cat2); 
        ind_cwcd = x1w==cat1(2)&x1s==cat1(2);
        ind_swsd = x1w==cat1(1)&x1s==cat1(1);
end

%% PCA
% run a PCA on ALL the data together, then plot the different subgroups of the data
% - base grouping: sediment core, sediment satellite, water core, water satellite
% - additional: plot MF that are both in core water & core sediment, core water & satellite sediment,
%                    core sediment & satellite water, satellite water & satellite sediment
X = [data.AI, data.AI_Mod, data.DBE, data.DBE_O, data.DBE_AI, data.GFE, data.NOSC, data.kdefectCH2, data.OtoC_ratio, data.HtoC_ratio, data.NtoC_ratio, data.PtoC_ratio, data.NtoC_ratio];
vbls={'AI','AI_Mod','DBE','DBE_O','DBE_AI','GFE','NOSC','kdefectCH2','O/C', 'H/C', 'N/C', 'P/C', 'N/C'};
Z = zscore(X); % Standardized data
[coefs,score] = pca(Z);
ndim = 2; % dimensions of the biplot

figure,
subplot(4,2,1), biplot(coefs(:,1:ndim),'Scores',score(x1w==cat1(1),1:ndim),'VarLabels',vbls); title(strcat(ti1,' water,satellite'))
subplot(4,2,3), biplot(coefs(:,1:ndim),'Scores',score(x1w==cat1(2),1:ndim),'VarLabels',vbls); title(strcat(ti1,' water,core'))
subplot(4,2,2), biplot(coefs(:,1:ndim),'Scores',score(x1s==cat1(1),1:ndim),'VarLabels',vbls); title(strcat(ti1,' sed,satellite'))
subplot(4,2,4), biplot(coefs(:,1:ndim),'Scores',score(x1s==cat1(2),1:ndim),'VarLabels',vbls); title(strcat(ti1,' sed,core'))
subplot(4,2,5), biplot(coefs(:,1:ndim),'Scores',score(ind_cwcd,1:ndim),'VarLabels',vbls); title(strcat(ti1,' core water & core sed'))
subplot(4,2,6), biplot(coefs(:,1:ndim),'Scores',score(ind_swsd,1:ndim),'VarLabels',vbls); title(strcat(ti1,' sat water & sat sed'))
if(method==0)
    subplot(4,2,7), biplot(coefs(:,1:ndim),'Scores',score(x1w==cat1(3),1:ndim),'VarLabels',vbls); title(strcat(ti1,' water,in-between'))
    subplot(4,2,8), biplot(coefs(:,1:ndim),'Scores',score(x1s==cat1(3),1:ndim),'VarLabels',vbls); title(strcat(ti1,' sed,in-between'))
end
