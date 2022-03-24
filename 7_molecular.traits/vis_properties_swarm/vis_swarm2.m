% visualize some of the aggregate properties for core and satellite species
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


%% visualize the data with violin or swarm plots - all data
switch method
    case 0      % emergent
        cat1 = ["In-between" "Satellite" "Core"];
        x1w = categorical(data.csflagemergent_water,cat1);  tiw = 'Water, emergent';
        x1s = categorical(data.csflagemergent_sed,cat1);    tis = 'Sediment, emergent';
    case 1      % pca
        cat1 = ["Satellite" "Core"];
        x1w = categorical(data.csflagemergent_water,cat1);  tiw = 'Water, PCA';
        x1s = categorical(data.csflagemergent_sed,cat1);    tis = 'Sediment, PCA';
    case 2      % random forest
        cat1 = ["Satellite" "Core"];
        x1w = categorical(data.csflagrf_water,cat1);        tiw = 'Water, random forest';
        x1s = categorical(data.csflagrf_sed,cat1);          tis = 'Sediment, random forest';
end

n=13;
figure1 = figure;
%  water
subplot(n,2, 1), swarmchart(x1w,data.DBE,'.'), title('DBE')
subplot(n,2, 3), swarmchart(x1w,data.DBE_O,'.'), title('DBE O')
subplot(n,2, 5), swarmchart(x1w,data.DBE_AI,'.'), title('DBE AI')
subplot(n,2, 7), swarmchart(x1w,data.AI,'.'), title('AI')
subplot(n,2, 9), swarmchart(x1w,data.AI_Mod,'.'), title('AI mod')
subplot(n,2,11), swarmchart(x1w,data.GFE,'.'), title('GFE')
subplot(n,2,13), swarmchart(x1w,data.NOSC,'.'), title('NOSC'), hold on
subplot(n,2,15), swarmchart(x1w,data.kdefectCH2,'.'), title('kdefectCH2'), hold on
subplot(n,2,17), swarmchart(x1w,data.OtoC_ratio,'.'), title('OtoC ratio'), hold on
subplot(n,2,19), swarmchart(x1w,data.HtoC_ratio,'.'), title('HtoC ratio'), hold on
subplot(n,2,21), swarmchart(x1w,data.NtoC_ratio,'.'), title('NtoC ratio'), hold on
subplot(n,2,23), swarmchart(x1w,data.PtoC_ratio,'.'), title('PtoC ratio'), hold on
subplot(n,2,25), swarmchart(x1w,data.NtoP_ratio,'.'), title('NtoP ratio'), hold on
annotation(figure1,'textbox', [0.13 0.94 0.16 0.02],'String',{tiw},'FitBoxToText','on');
%  sediment
subplot(n,2, 2), swarmchart(x1s,data.DBE,'.'), title('DBE')
subplot(n,2, 4), swarmchart(x1s,data.DBE_O,'.'), title('DBE O')
subplot(n,2, 6), swarmchart(x1s,data.DBE_AI,'.'), title('DBE AI')
subplot(n,2, 8), swarmchart(x1s,data.AI,'.'), title('AI')
subplot(n,2,10), swarmchart(x1s,data.AI_Mod,'.'), title('AI mod')
subplot(n,2,12), swarmchart(x1s,data.GFE,'.'), title('GFE')
subplot(n,2,14), swarmchart(x1s,data.NOSC,'.'), title('NOSC'), hold on
subplot(n,2,16), swarmchart(x1s,data.kdefectCH2,'.'), title('kdefectCH2'), hold on
subplot(n,2,18), swarmchart(x1s,data.OtoC_ratio,'.'), title('OtoC ratio'), hold on
subplot(n,2,20), swarmchart(x1s,data.HtoC_ratio,'.'), title('HtoC ratio'), hold on
subplot(n,2,22), swarmchart(x1s,data.NtoC_ratio,'.'), title('NtoC ratio'), hold on
subplot(n,2,24), swarmchart(x1s,data.PtoC_ratio,'.'), title('PtoC ratio'), hold on
subplot(n,2,26), swarmchart(x1s,data.NtoP_ratio,'.'), title('NtoP ratio'), hold on
annotation(figure1,'textbox', [0.57 0.94 0.16 0.02],'String',{tis},'FitBoxToText','on');

%% visualize the data after grouping them into compound classes

catclass=unique(data.Class);
xclass = categorical(data.Class,catclass);

% loop over compound classes
for icat = 1:size(catclass,1)

    % select MF in a given compound class
    ind = (data.Class==catclass(icat));
    data2= data(ind,:);

    switch method
        case 0      % emergent
            cat1 = ["In-between" "Satellite" "Core"];
            x1w = categorical(data2.csflagemergent_water,cat1);  tiw = 'Water, emergent';
            x1s = categorical(data2.csflagemergent_sed,cat1);    tis = 'Sediment, emergent';
        case 1      % pca
            cat1 = ["Satellite" "Core"];
            x1w = categorical(data2.csflagemergent_water,cat1);  tiw = 'Water, PCA';
            x1s = categorical(data2.csflagemergent_sed,cat1);    tis = 'Sediment, PCA';
        case 2      % random forest
            cat1 = ["Satellite" "Core"];
            x1w = categorical(data2.csflagrf_water,cat1);        tiw = 'Water, random forest';
            x1s = categorical(data2.csflagrf_sed,cat1);          tis = 'Sediment, random forest';
    end

    n=13;
    figure1 = figure;
    %  water
    subplot(n,2, 1), swarmchart(x1w,data2.DBE,'.'), title('DBE')
    subplot(n,2, 3), swarmchart(x1w,data2.DBE_O,'.'), title('DBE O')
    subplot(n,2, 5), swarmchart(x1w,data2.DBE_AI,'.'), title('DBE AI')
    subplot(n,2, 7), swarmchart(x1w,data2.AI,'.'), title('AI')
    subplot(n,2, 9), swarmchart(x1w,data2.AI_Mod,'.'), title('AI mod')
    subplot(n,2,11), swarmchart(x1w,data2.GFE,'.'), title('GFE')
    subplot(n,2,13), swarmchart(x1w,data2.NOSC,'.'), title('NOSC'), hold on
    subplot(n,2,15), swarmchart(x1w,data2.kdefectCH2,'.'), title('kdefectCH2'), hold on
    subplot(n,2,17), swarmchart(x1w,data2.OtoC_ratio,'.'), title('OtoC ratio'), hold on
    subplot(n,2,19), swarmchart(x1w,data2.HtoC_ratio,'.'), title('HtoC ratio'), hold on
    subplot(n,2,21), swarmchart(x1w,data2.NtoC_ratio,'.'), title('NtoC ratio'), hold on
    subplot(n,2,23), swarmchart(x1w,data2.PtoC_ratio,'.'), title('PtoC ratio'), hold on
    subplot(n,2,25), swarmchart(x1w,data2.NtoP_ratio,'.'), title('NtoP ratio'), hold on
    annotation(figure1,'textbox', [0.13 0.94 0.16 0.02],'String',{strcat(tiw, string(catclass(icat)))},'FitBoxToText','on');
    %  sediment
    subplot(n,2, 2), swarmchart(x1s,data2.DBE,'.'), title('DBE')
    subplot(n,2, 4), swarmchart(x1s,data2.DBE_O,'.'), title('DBE O')
    subplot(n,2, 6), swarmchart(x1s,data2.DBE_AI,'.'), title('DBE AI')
    subplot(n,2, 8), swarmchart(x1s,data2.AI,'.'), title('AI')
    subplot(n,2,10), swarmchart(x1s,data2.AI_Mod,'.'), title('AI mod')
    subplot(n,2,12), swarmchart(x1s,data2.GFE,'.'), title('GFE')
    subplot(n,2,14), swarmchart(x1s,data2.NOSC,'.'), title('NOSC'), hold on
    subplot(n,2,16), swarmchart(x1s,data2.kdefectCH2,'.'), title('kdefectCH2'), hold on
    subplot(n,2,18), swarmchart(x1s,data2.OtoC_ratio,'.'), title('OtoC ratio'), hold on
    subplot(n,2,20), swarmchart(x1s,data2.HtoC_ratio,'.'), title('HtoC ratio'), hold on
    subplot(n,2,22), swarmchart(x1s,data2.NtoC_ratio,'.'), title('NtoC ratio'), hold on
    subplot(n,2,24), swarmchart(x1s,data2.PtoC_ratio,'.'), title('PtoC ratio'), hold on
    subplot(n,2,26), swarmchart(x1s,data2.NtoP_ratio,'.'), title('NtoP ratio'), hold on
    annotation(figure1,'textbox', [0.57 0.94 0.16 0.02],'String',{strcat(tis, string(catclass(icat)))},'FitBoxToText','on');

end


%% visualize core water&sed, satellite water&sed
switch method
    case 0      % emergent
        cat1 = ["Global core" "Global satellite","Shifter"];
        x1 = categorical(data.csflagemergent_generaloverlap,cat1);  ti = 'global core or satellite: emergent';
    case 1      % pca
        cat1 = ["Global core" "Global satellite","Shifter"];
        x1 = categorical(data.csflagpca_generaloverlap,cat1);  ti = 'global core or satellite: pca';
    case 2      % random forest
        cat1 = ["Global core" "Global satellite","Shifter"];
        x1 = categorical(data.csflagrf_generaloverlap,cat1);  ti = 'global core or satellite: random forest';
end

n=13;
figure1 = figure;
subplot(n,1, 1), swarmchart(x1,data.DBE,'.'), title('DBE')
subplot(n,1, 2), swarmchart(x1,data.DBE_O,'.'), title('DBE O')
subplot(n,1, 3), swarmchart(x1,data.DBE_AI,'.'), title('DBE AI')
subplot(n,1, 4), swarmchart(x1,data.AI,'.'), title('AI')
subplot(n,1, 5), swarmchart(x1,data.AI_Mod,'.'), title('AI mod')
subplot(n,1, 6), swarmchart(x1,data.GFE,'.'), title('GFE')
subplot(n,1, 7), swarmchart(x1,data.NOSC,'.'), title('NOSC'), hold on
subplot(n,1, 8), swarmchart(x1,data.kdefectCH2,'.'), title('kdefectCH2'), hold on
subplot(n,1, 9), swarmchart(x1,data.OtoC_ratio,'.'), title('OtoC ratio'), hold on
subplot(n,1,10), swarmchart(x1,data.HtoC_ratio,'.'), title('HtoC ratio'), hold on
subplot(n,1,11), swarmchart(x1,data.NtoC_ratio,'.'), title('NtoC ratio'), hold on
subplot(n,1,12), swarmchart(x1,data.PtoC_ratio,'.'), title('PtoC ratio'), hold on
subplot(n,1,13), swarmchart(x1,data.NtoP_ratio,'.'), title('NtoP ratio'), hold on
annotation(figure1,'textbox', [0.13 0.94 0.16 0.02],'String',{ti},'FitBoxToText','on');
