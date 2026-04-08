
% ADD 'doseResponse_RH.m' function to the path

clc
clear
warning off
close all

% EXCEL FILE DIRECTORY
directory = '/Users/Pablo/Documents/Bellono Lab/Data/InformationCoding/v2' ;
cd(directory)

%% LIGAND LIST TO FIND DATA IN EXCEL FILE
% Responses are normalized to the vehicle (no ligand)
ligandList = {
    'Androsterone'
    'Progesterone'
    'Chloroquine'
    'Strychnine'
    'Naringin'
    'Nootkatone'
    'Taurocholic acid'
    
    } ;


% Saturation dose for each molecule (in 'ligandList' order)
% Above saturation the response of the receptor goes down
maxDose_CR518     = [ 11 11 11 10 11  9  11  ]  ; % Blue
maxDose_CR918     = [ 11 11 11 10 11 11  11  ]  ; % Green
maxDose_CR999     = [ 11 11 11 10 11  8  11  ]  ; % Cyan
maxDose_CR518_918 = [ 11 11 11 10 11 11  11  ]  ; % Red
maxDose_CR518_999 = [ 11 11 11  9  8  8  11  ]  ; % Magenta

% ABOVE THIS THRESHOLD WILL ATTEMPT FITTING A SIGMOID
respThresh = 4 ; % arbitrary


%% GET ALL VALUES FROM EXCEL FILE
tableRaw_CR518     = readtable('doseHeteromer_dataAll_v2.xlsx', 'Sheet','CR518', 'VariableNamingRule','preserve', 'Range','A1:Z1000')      ;
tableRaw_CR918     = readtable('doseHeteromer_dataAll_v2.xlsx', 'Sheet','CR918', 'VariableNamingRule','preserve', 'Range','A1:Z1000')      ;
tableRaw_CR999     = readtable('doseHeteromer_dataAll_v2.xlsx', 'Sheet','CR999', 'VariableNamingRule','preserve', 'Range','A1:Z1000')      ;
tableRaw_CR518_918 = readtable('doseHeteromer_dataAll_v2.xlsx', 'Sheet','CR518-918', 'VariableNamingRule','preserve', 'Range','A1:Z1000')  ;
tableRaw_CR518_999 = readtable('doseHeteromer_dataAll_v2.xlsx', 'Sheet','CR518-999', 'VariableNamingRule','preserve', 'Range','A1:Z1000')  ;

% GCaMP: Cells express only the Ca-sensor to monitor endogenous responses that are independent of CR expression
tableRaw_GCaMP     = readtable('doseHeteromer_dataAll_v2.xlsx', 'Sheet','GCaMP', 'VariableNamingRule','preserve', 'Range','A1:Z1000')      ;

agonistIdx_CR518 = [ ] ;
for nAgonist= 1: length(ligandList)
    agonistIdx_CR518(:, nAgonist)   =  find(strcmp(tableRaw_CR518.Var2, ligandList(nAgonist)))  ;
end

agonistIdx_CR918 = [ ] ;
for nAgonist= 1: length(ligandList)
    agonistIdx_CR918(:, nAgonist)   =  find(strcmp(tableRaw_CR918.Var2, ligandList(nAgonist)))  ;
end

agonistIdx_CR999 = [ ] ;
for nAgonist= 1: length(ligandList)
    agonistIdx_CR999(:, nAgonist)   =  find(strcmp(tableRaw_CR999.Var2, ligandList(nAgonist)))  ;
end

agonistIdx_CR518_918 = [ ] ;
for nAgonist= 1: length(ligandList)
    agonistIdx_CR518_918(:, nAgonist)   =  find(strcmp(tableRaw_CR518_918.Var2, ligandList(nAgonist)))  ;
end

agonistIdx_CR518_999 = [ ] ;
for nAgonist= 1: length(ligandList)
    agonistIdx_CR518_999(:, nAgonist)   =  find(strcmp(tableRaw_CR518_999.Var2, ligandList(nAgonist)))  ;
end

agonistIdx_GCaMP = [ ] ;
for nAgonist= 1: length(ligandList)
    agonistIdx_GCaMP(:, nAgonist)   =  find(strcmp(tableRaw_GCaMP.Var2, ligandList(nAgonist)))  ;
end



%% GET VALUES FOR EACH RECEPTOR
CR518_all= [] ;
for nAgonist= 1: length(ligandList)
    tempValues = table2array(tableRaw_CR518(agonistIdx_CR518(:, nAgonist), 4:end)) ; % 3=first, 12=last concentration
    CR518_all{nAgonist, 1} = tempValues ;
end

CR918_all= [] ;
for nAgonist= 1: length(ligandList)
    tempValues = table2array(tableRaw_CR918(agonistIdx_CR918(:, nAgonist), 4:end)) ; % 3=first, 12=last concentration
    CR918_all{nAgonist, 1} = tempValues ;
end

CR999_all= [] ;
for nAgonist= 1: length(ligandList)
    tempValues = table2array(tableRaw_CR999(agonistIdx_CR999(:, nAgonist), 4:end)) ; % 3=first, 12=last concentration
    CR999_all{nAgonist, 1} = tempValues ;
end

CR518_918_all= [] ;
for nAgonist= 1: length(ligandList)
    tempValues = table2array(tableRaw_CR518_918(agonistIdx_CR518_918(:, nAgonist), 4:end)) ; % 3=first, 12=last concentration
    CR518_918_all{nAgonist, 1} = tempValues ;
end

CR518_999_all= [] ;
for nAgonist= 1: length(ligandList)
    tempValues = table2array(tableRaw_CR518_999(agonistIdx_CR518_999(:, nAgonist), 4:end)) ; % 3=first, 12=last concentration
    CR518_999_all{nAgonist, 1} = tempValues ;
end

GCaMP_all= [] ;
for nAgonist= 1: length(ligandList)
    tempValues = table2array(tableRaw_GCaMP(agonistIdx_GCaMP(:, nAgonist), 4:end)) ; % 3=first, 12=last concentration
    GCaMP_all{nAgonist, 1} = tempValues ;
end

% PLOT ALL RECEPTORS/SINGLE MOLECULE
figure
tiling= tiledlayout(3,5, "TileSpacing","compact") ;
xlabel(tiling,'[Agonist] (µM)', 'FontSize', 16)
ylabel(tiling,'Vehicle normalized response', 'FontSize', 16)

EC50 = [ ] ;
for nAgonist= 1: length(ligandList)
    nexttile
    CR518_Mean_temp     = mean(CR518_all{nAgonist, 1}, 2, 'omitnan')     ;
    CR518_918_Mean_temp = mean(CR518_918_all{nAgonist, 1}, 2, 'omitnan') ;
    CR518_999_Mean_temp = mean(CR518_999_all{nAgonist, 1}, 2, 'omitnan') ;
    GCaMP_Mean_temp     = mean(GCaMP_all{nAgonist, 1}, 2, 'omitnan')     ;
    CR918_Mean_temp     = mean(CR918_all{nAgonist, 1}, 2, 'omitnan')     ;
    CR999_Mean_temp     = mean(CR999_all{nAgonist, 1}, 2, 'omitnan')     ;

    CR518_SEM_temp      = std(CR518_all{nAgonist, 1}, 0, 2, 'omitnan')./sqrt(sum(~isnan(CR518_all{nAgonist, 1}), 2))         ;
    CR518_918_SEM_temp  = std(CR518_918_all{nAgonist, 1}, 0, 2, 'omitnan')./sqrt(sum(~isnan(CR518_918_all{nAgonist, 1}), 2)) ;
    CR518_999_SEM_temp  = std(CR518_999_all{nAgonist, 1}, 0, 2, 'omitnan')./sqrt(sum(~isnan(CR518_999_all{nAgonist, 1}), 2)) ;
    GCaMP_SEM_temp      = std(GCaMP_all{nAgonist, 1}, 0, 2, 'omitnan')./sqrt(sum(~isnan(GCaMP_all{nAgonist, 1}), 2))         ;
    CR918_SEM_temp      = std(CR918_all{nAgonist, 1}, 0, 2, 'omitnan')./sqrt(sum(~isnan(CR918_all{nAgonist, 1}), 2))         ;
    CR999_SEM_temp      = std(CR999_all{nAgonist, 1}, 0, 2, 'omitnan')./sqrt(sum(~isnan(CR999_all{nAgonist, 1}), 2))         ;

    % CR518 BLUE
    doseTempCR518 = table2array(tableRaw_CR518(agonistIdx_CR518(:, nAgonist), 'Var3'));

    if any(CR518_Mean_temp > respThresh) == 1
        [hillCoefftmp, ec50tmp_CR518, minDose, maxDoseOut, coeffs, meanResponse, doses, sigmoid] = ...
            doseResponse_RH(doseTempCR518(1:maxDose_CR518(nAgonist)),CR518_Mean_temp(1:maxDose_CR518(nAgonist)), 1) ;

        EC50(nAgonist).CR518      = real(ec50tmp_CR518) ;
        hillCoeff = hillCoefftmp  ;
        xpoints   =logspace(log10(minDose),log10(maxDoseOut),1000) ;
        semilogx(xpoints,sigmoid(coeffs,xpoints),'Color','b','LineWidth', 2) ;
        hold on
        semilogx(doseTempCR518(1:maxDose_CR518(nAgonist)),CR518_Mean_temp(1:maxDose_CR518(nAgonist)), 'LineStyle','none') ;
        hold on
        errorbar(doseTempCR518(1:maxDose_CR518(nAgonist)),CR518_Mean_temp(1:maxDose_CR518(nAgonist)), [], CR518_SEM_temp(1:maxDose_CR518(nAgonist)),...
            'o','Color','k','MarkerSize', 5,'MarkerFaceColor','b','LineStyle','none', 'LineWidth', 1) ;

        text(0.05, 0.9, [sprintf('EC_{50}=%0.3g', EC50(nAgonist).CR518) ' µM'], ...
            'Units', 'normalized', 'FontSize', 10, 'FontName', 'Arial', 'FontWeight', 'bold', ...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Color','b');
    else
        errorbar(doseTempCR518(1:maxDose_CR518(nAgonist)),CR518_Mean_temp(1:maxDose_CR518(nAgonist)), [], CR518_SEM_temp(1:maxDose_CR518(nAgonist)),...
            'o','Color','k','MarkerSize', 5,'MarkerFaceColor','b','LineStyle','none', 'LineWidth', 1) ;

        text(0.05, 0.9, 'EC_{50} = [ ]', ...
            'Units', 'normalized', 'FontSize', 10, 'FontName', 'Arial', 'FontWeight', 'bold', ...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Color','b');
    end

    % CR518-918 RED
    doseTemp_518_918 = table2array(tableRaw_CR518_918(agonistIdx_CR518_918(:, nAgonist), 'Var3'));

    if any(CR518_918_Mean_temp > respThresh) == 1
        [hillCoefftmp, ec50tmp_CR518_918, minDose, maxDoseOut, coeffs, meanResponse, doses, sigmoid] = ...
            doseResponse_RH(doseTemp_518_918(1:maxDose_CR518_918(nAgonist)),CR518_918_Mean_temp(1:maxDose_CR518_918(nAgonist)), 1) ;

        EC50(nAgonist).CR518_918      = real(ec50tmp_CR518_918) ;
        hillCoeff = hillCoefftmp  ;
        xpoints   =logspace(log10(minDose),log10(maxDoseOut),1000) ;
        semilogx(xpoints,sigmoid(coeffs,xpoints),'Color','r','LineWidth', 2) ;
        hold on
        semilogx(doseTemp_518_918(1:maxDose_CR518_918(nAgonist)),CR518_918_Mean_temp(1:maxDose_CR518_918(nAgonist)), 'LineStyle','none') ;
        hold on
        errorbar(doseTemp_518_918(1:maxDose_CR518_918(nAgonist)),CR518_918_Mean_temp(1:maxDose_CR518_918(nAgonist)), [], CR518_918_SEM_temp(1:maxDose_CR518_918(nAgonist)),...
            'o','Color','k','MarkerSize', 5,'MarkerFaceColor','r','LineStyle','none', 'LineWidth', 1) ;

        text(0.05, 0.8, [sprintf('EC_{50}=%0.3g', EC50(nAgonist).CR518_918 ) ' µM'], ...
            'Units', 'normalized', 'FontSize', 10, 'FontName', 'Arial', 'FontWeight', 'bold', ...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Color','r');
    else
        errorbar(doseTemp_518_918(1:maxDose_CR518_918(nAgonist)),CR518_918_Mean_temp(1:maxDose_CR518_918(nAgonist)), [], CR518_918_SEM_temp(1:maxDose_CR518_918(nAgonist)),...
            'o','Color','k','MarkerSize', 5,'MarkerFaceColor','r','LineStyle','none', 'LineWidth', 1) ;

        text(0.05, 0.8, 'EC_{50} = [ ]', ...
            'Units', 'normalized', 'FontSize', 10, 'FontName', 'Arial', 'FontWeight', 'bold', ...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Color','r');
    end

    % CR918 GREEN
    doseTempCR918 = table2array(tableRaw_CR918(agonistIdx_CR918(:, nAgonist), 'Var3'));

    if any(CR918_Mean_temp > respThresh) == 1
        [hillCoefftmp, ec50tmp_CR918, minDose, maxDoseOut, coeffs, meanResponse, doses, sigmoid] = ...
            doseResponse_RH(doseTempCR918(1:maxDose_CR918(nAgonist)),CR918_Mean_temp(1:maxDose_CR918(nAgonist)), 1) ;

        EC50(nAgonist).CR918      = real(ec50tmp_CR918) ;
        hillCoeff = hillCoefftmp  ;
        xpoints   =logspace(log10(minDose),log10(maxDoseOut),1000) ;
        semilogx(xpoints,sigmoid(coeffs,xpoints),'Color','g','LineWidth', 2) ;
        hold on
        semilogx(doseTempCR918(1:maxDose_CR918(nAgonist)),CR918_Mean_temp(1:maxDose_CR918(nAgonist)), 'LineStyle','none') ;
        hold on
        errorbar(doseTempCR918(1:maxDose_CR918(nAgonist)),CR918_Mean_temp(1:maxDose_CR918(nAgonist)), [], CR918_SEM_temp(1:maxDose_CR918(nAgonist)),...
            'o','Color','k','MarkerSize', 5,'MarkerFaceColor','g','LineStyle','none', 'LineWidth', 1) ;

        text(0.05, 0.7, [sprintf('EC_{50}=%0.3g', EC50(nAgonist).CR918 ) ' µM'], ...
            'Units', 'normalized', 'FontSize', 10, 'FontName', 'Arial', 'FontWeight', 'bold', ...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Color','g');
    else
        errorbar(doseTempCR918(1:maxDose_CR918(nAgonist)),CR918_Mean_temp(1:maxDose_CR918(nAgonist)), [], CR918_SEM_temp(1:maxDose_CR918(nAgonist)),...
            'o','Color','k','MarkerSize', 5,'MarkerFaceColor','g','LineStyle','none', 'LineWidth', 1) ;

        text(0.05, 0.7, 'EC_{50} = [ ]', ...
            'Units', 'normalized', 'FontSize', 10, 'FontName', 'Arial', 'FontWeight', 'bold', ...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Color','g');
    end

    % CR999 CYAN
    doseTempCR999 = table2array(tableRaw_CR999(agonistIdx_CR999(:, nAgonist), 'Var3'));

    if any(CR999_Mean_temp > respThresh) == 1
        [hillCoefftmp, ec50tmp_CR999, minDose, maxDoseOut, coeffs, meanResponse, doses, sigmoid] = ...
            doseResponse_RH(doseTempCR999(1:maxDose_CR999(nAgonist)),CR999_Mean_temp(1:maxDose_CR999(nAgonist)), 1) ;

        EC50(nAgonist).CR999      = real(ec50tmp_CR999) ;
        hillCoeff = hillCoefftmp  ;
        xpoints   =logspace(log10(minDose),log10(maxDoseOut),1000) ;
        semilogx(xpoints,sigmoid(coeffs,xpoints),'Color','c','LineWidth', 2) ;
        hold on
        semilogx(doseTempCR999(1:maxDose_CR999(nAgonist)),CR999_Mean_temp(1:maxDose_CR999(nAgonist)), 'LineStyle','none') ;
        hold on
        errorbar(doseTempCR999(1:maxDose_CR999(nAgonist)),CR999_Mean_temp(1:maxDose_CR999(nAgonist)), [], CR999_SEM_temp(1:maxDose_CR999(nAgonist)),...
            'o','Color','k','MarkerSize', 5,'MarkerFaceColor','c','LineStyle','none', 'LineWidth', 1) ;

        text(0.05, 0.6, [sprintf('EC_{50}=%0.3g', EC50(nAgonist).CR999 ) ' µM'], ...
            'Units', 'normalized', 'FontSize', 10, 'FontName', 'Arial', 'FontWeight', 'bold', ...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Color','c');
    else
        errorbar(doseTempCR999(1:maxDose_CR999(nAgonist)),CR999_Mean_temp(1:maxDose_CR999(nAgonist)), [], CR999_SEM_temp(1:maxDose_CR999(nAgonist)),...
            'o','Color','k','MarkerSize', 5,'MarkerFaceColor','c','LineStyle','none', 'LineWidth', 1) ;

        text(0.05, 0.6, 'EC_{50} = [ ]', ...
            'Units', 'normalized', 'FontSize', 10, 'FontName', 'Arial', 'FontWeight', 'bold', ...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Color','c');
    end

    % CR518-999 MAGENTA
    doseTemp_518_999 = table2array(tableRaw_CR518_999(agonistIdx_CR518_999(:, nAgonist), 'Var3'));

    if any(CR518_999_Mean_temp > respThresh) == 1
        [hillCoefftmp, ec50tmp_CR518_999, minDose, maxDoseOut, coeffs, meanResponse, doses, sigmoid] = ...
            doseResponse_RH(doseTemp_518_999(1:maxDose_CR518_999(nAgonist)),CR518_999_Mean_temp(1:maxDose_CR518_999(nAgonist)), 1) ;

        EC50(nAgonist).CR518_999      = real(ec50tmp_CR518_999) ;
        hillCoeff = hillCoefftmp  ;
        xpoints   =logspace(log10(minDose),log10(maxDoseOut),1000) ;
        semilogx(xpoints,sigmoid(coeffs,xpoints),'Color','m','LineWidth', 2) ;
        hold on
        semilogx(doseTemp_518_999(1:maxDose_CR518_999(nAgonist)),CR518_999_Mean_temp(1:maxDose_CR518_999(nAgonist)), 'LineStyle','none') ;
        hold on
        errorbar(doseTemp_518_999(1:maxDose_CR518_999(nAgonist)),CR518_999_Mean_temp(1:maxDose_CR518_999(nAgonist)), [], CR518_999_SEM_temp(1:maxDose_CR518_999(nAgonist)),...
            'o','Color','k','MarkerSize', 5,'MarkerFaceColor','m','LineStyle','none', 'LineWidth', 1) ;

        text(0.05, 0.5, [sprintf('EC_{50}=%0.3g', EC50(nAgonist).CR518_999 ) ' µM'], ...
            'Units', 'normalized', 'FontSize', 10, 'FontName', 'Arial', 'FontWeight', 'bold', ...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Color','m');
    else
        errorbar(doseTemp_518_999(1:maxDose_CR518_999(nAgonist)),CR518_999_Mean_temp(1:maxDose_CR518_999(nAgonist)), [], CR518_999_SEM_temp(1:maxDose_CR518_999(nAgonist)),...
            'o','Color','k','MarkerSize', 5,'MarkerFaceColor','m','LineStyle','none', 'LineWidth', 1) ;

        text(0.05, 0.5, 'EC_{50} = [ ]', ...
            'Units', 'normalized', 'FontSize', 10, 'FontName', 'Arial', 'FontWeight', 'bold', ...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Color','m');
    end

    % Find axes limits across receptors for each ligand/plot:
    maxDose_allReceptors = max([doseTempCR518(maxDose_CR518(nAgonist))
        doseTemp_518_918(maxDose_CR518_918(nAgonist))
        doseTempCR999(maxDose_CR999(nAgonist))
        doseTemp_518_999(maxDose_CR518_999(nAgonist))
        doseTempCR918(maxDose_CR918(nAgonist))]) ;

    minDose_allReceptors = min([doseTempCR518(1)
        doseTemp_518_918(1)
        doseTempCR999(1)
        doseTemp_518_999(1)
        doseTempCR918(1)])  ;

    if isnan(max(CR518_SEM_temp))== 1
        yMax_margin = 5 ;
    else
        yMax_margin= max(CR518_SEM_temp) ;
    end

    % yMax= round(max(max([CR518_Mean_temp(1:maxDose_CR518(nAgonist))...
    %                      CR518_918_Mean_temp(1:maxDose_CR518_918(nAgonist))...
    %                      CR918_Mean_temp(1:maxDose_CR918(nAgonist))]))) + yMax_margin ;


    yMax= round(max([max(CR518_Mean_temp(1:maxDose_CR518(nAgonist)))...
        max(CR518_918_Mean_temp(1:maxDose_CR518_918(nAgonist)))...
        max(CR999_Mean_temp(1:maxDose_CR999(nAgonist)))...
        max(CR518_999_Mean_temp(1:maxDose_CR518_999(nAgonist)))...
        max(CR918_Mean_temp(1:maxDose_CR918(nAgonist))) ])) + yMax_margin ;

    ylim([0 yMax + 1]), xlim([minDose_allReceptors maxDose_allReceptors]), title(ligandList(nAgonist)), box off
    set(gca, 'FontSize', 12, 'FontName', 'Arial');
    EC50(nAgonist).ligands= ligandList{nAgonist} ;

    set(gca, 'XScale', 'log')

end

% add legend
h1 = plot(nan, nan, 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k', 'DisplayName', 'CR518');
h2 = plot(nan, nan, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'DisplayName', 'CR518-918');
h3 = plot(nan, nan, 'o', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k', 'DisplayName', 'CR918');
h4 = plot(nan, nan, 'o', 'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'k', 'DisplayName', 'CR999');
h5 = plot(nan, nan, 'o', 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'k', 'DisplayName', 'CR518-999');

lgd = legend([h1 h2 h3 h4 h5], 'Orientation', 'horizontal', 'FontSize', 14);
lgd.Layout.Tile = 'south';  % or 'north', 'east', 'west'

set(gcf, 'Position',  [500, 500, 1500, 800]) % % [pos, pos, width, height]




%% MAX AMPLITUD HEATMAP
receptors = {'CR518', 'CR518-918', 'CR999', 'CR518-999', 'CR918'} ;

% Use maxDose for each agonist to extract values at max concentration
stackMeans = zeros(length(ligandList), 5) ;
stackSEMs  = zeros(length(ligandList), 5) ;

for nligand = 1:length(ligandList)
    % idx = maxDose_CR518(nligand);
    stackMeans(nligand,1) = mean(CR518_all{nligand}(maxDose_CR518(nligand),:), 'omitnan');
    stackMeans(nligand,2) = mean(CR518_918_all{nligand}(maxDose_CR518_918(nligand),:), 'omitnan');
    stackMeans(nligand,3) = mean(CR999_all{nligand}(maxDose_CR999(nligand),:), 'omitnan');
    stackMeans(nligand,4) = mean(CR518_999_all{nligand}(maxDose_CR518_999(nligand),:), 'omitnan');
    stackMeans(nligand,5) = mean(CR918_all{nligand}(maxDose_CR918(nligand),:), 'omitnan');

    stackSEMs(nligand,1) = std(CR518_all{nligand}(maxDose_CR518(nligand),:), 0, 2, 'omitnan') ./ sqrt(sum(~isnan(CR518_all{nligand}(maxDose_CR518(nligand),:))))   ;
    stackSEMs(nligand,2) = std(CR518_918_all{nligand}(maxDose_CR518_918(nligand),:), 0, 2, 'omitnan') ./ sqrt(sum(~isnan(CR518_918_all{nligand}(maxDose_CR518_918(nligand),:))))    ;
    stackSEMs(nligand,3) = std(CR999_all{nligand}(maxDose_CR999(nligand),:), 0, 2, 'omitnan') ./ sqrt(sum(~isnan(CR999_all{nligand}(maxDose_CR999(nligand),:))))             ;
    stackSEMs(nligand,4) = std(CR518_999_all{nligand}(maxDose_CR518_999(nligand),:), 0, 2, 'omitnan') ./ sqrt(sum(~isnan(CR518_999_all{nligand}(maxDose_CR518_999(nligand),:))))             ;
    stackSEMs(nligand,5) = std(CR918_all{nligand}(maxDose_CR918(nligand),:), 0, 2, 'omitnan') ./ sqrt(sum(~isnan(CR918_all{nligand}(maxDose_CR918(nligand),:))))             ;
end

figure
h = heatmap(ligandList, receptors, stackMeans', ...
    'Colormap', parula, ...
    'MissingDataColor', 'w', ...
    'MissingDataLabel', 'N/A');

title('Peak response (vehicle folds)')
h.CellLabelColor = 'none';  % Hide numbers
set(gca, 'FontSize', 12, 'FontName', 'Arial');
set(gcf, 'Position',  [500, 500, 350, 200]) % % [pos, pos, width, height]


%% EC50 HEATMAP
% Updated EC50 values (CR518 in Column1, CR518_918 in Column2)
% ENTER VALUES MANUALLY

EC50_valuesAll = nan(length(ligandList), 3);
for nVal = 1:length(ligandList)
    if ~isempty(EC50(nVal).CR518)
        EC50_valuesAll(nVal,1) = EC50(nVal).CR518;
    end
    if ~isempty(EC50(nVal).CR518_918)
        EC50_valuesAll(nVal,2) = EC50(nVal).CR518_918;
    end
    if ~isempty(EC50(nVal).CR918)
        EC50_valuesAll(nVal,3) = EC50(nVal).CR918;
    end
    if ~isempty(EC50(nVal).CR999)
        EC50_valuesAll(nVal,4) = EC50(nVal).CR999;
    end
    if ~isempty(EC50(nVal).CR518_999)
        EC50_valuesAll(nVal,5) = EC50(nVal).CR518_999;
    end
end

% Log-transform EC50
log_ec50  = log10(EC50_valuesAll)   ;
log_ec50  = real(log_ec50)          ;

figure
h = heatmap(ligandList, receptors, log_ec50', ...
    'Colormap', parula, ...
    'MissingDataColor', 'w', ...
    'MissingDataLabel', 'N/A');

title('log_{10}(EC50)')
h.CellLabelColor = 'none';  % Hide numbers
set(gca, 'FontSize', 12, 'FontName', 'Arial');
set(gcf, 'Position',  [500, 500, 350, 200]) % % [pos, pos, width, height]







