
% ADD 'doseResponse_RH.m' function to the path

clc
clear
warning off
close all

% EXCEL FILE DIRECTORY
directory = '/Users/Pablo/Documents/Bellono Lab/Data/Screenings/doseCurves/Heteromers' ;
cd(directory)

%% LIGAND LIST TO FIND DATA IN EXCEL FILE
% Responses are normalized to the vehicle (no ligand)
ligandList = {
    'Medroxyprogesterone'
    'Androsterone'
    'Progesterone'
    'Estrone'
    'Chloroquine'
    'Strychnine'
    'Naringin'
    'Nootkatone'
    'Costunolide'
    'Taurocholic acid'
    'Norharmane' } ;

% Saturation dose for each molecule (in 'ligandList' order)
% Above saturation the response of the receptor goes down
maxDose_CR518     = [9 11 11 11 11 10 11 9 8  11 11  ] ;
maxDose_CR518_918 = [9 11 11 11 11 10 11 11 11 10 11  ] ;
maxDose_CR918     = [9 11 11 11 11 10 11 11 11 10 11  ] ;


% ABOVE THIS THRESHOLD WILL ATTEMPT FITTING A SIGMOID 
respThresh = 3 ;


%% GET ALL VALUES FROM EXCEL FILE
tableRaw_CR518     = readtable('doseHeteromer_dataAll.xlsx', 'Sheet','CR518', 'VariableNamingRule','preserve', 'Range','A1:Z1000')      ;
tableRaw_CR518_918 = readtable('doseHeteromer_dataAll.xlsx', 'Sheet','CR518-918', 'VariableNamingRule','preserve', 'Range','A1:Z1000')  ;
tableRaw_CR918     = readtable('doseHeteromer_dataAll.xlsx', 'Sheet','CR918', 'VariableNamingRule','preserve', 'Range','A1:Z1000')      ;

% GCaMP: Cells express only the Ca-sensor to monitor endogenous responses that are independent of CR expression
tableRaw_GCaMP     = readtable('doseHeteromer_dataAll.xlsx', 'Sheet','GCaMP', 'VariableNamingRule','preserve', 'Range','A1:Z1000')      ;

agonistIdx_CR518 = [ ] ;
for nAgonist= 1: length(ligandList)
    agonistIdx_CR518(:, nAgonist)   =  find(strcmp(tableRaw_CR518.Var2, ligandList(nAgonist)))  ;
end

agonistIdx_CR518_918 = [ ] ;
for nAgonist= 1: length(ligandList)
    agonistIdx_CR518_918(:, nAgonist)   =  find(strcmp(tableRaw_CR518_918.Var2, ligandList(nAgonist)))  ;
end

agonistIdx_CR918 = [ ] ;
for nAgonist= 1: length(ligandList)
    agonistIdx_CR918(:, nAgonist)   =  find(strcmp(tableRaw_CR918.Var2, ligandList(nAgonist)))  ;
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

CR518_918_all= [] ;
for nAgonist= 1: length(ligandList)
    tempValues = table2array(tableRaw_CR518_918(agonistIdx_CR518_918(:, nAgonist), 4:end)) ; % 3=first, 12=last concentration
    CR518_918_all{nAgonist, 1} = tempValues ;
end

CR918_all= [] ;
for nAgonist= 1: length(ligandList)
    tempValues = table2array(tableRaw_CR918(agonistIdx_CR518_918(:, nAgonist), 4:end)) ; % 3=first, 12=last concentration
    CR918_all{nAgonist, 1} = tempValues ;
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
    GCaMP_Mean_temp     = mean(GCaMP_all{nAgonist, 1}, 2, 'omitnan')     ;
    CR918_Mean_temp     = mean(CR918_all{nAgonist, 1}, 2, 'omitnan')     ;

    CR518_SEM_temp      = std(CR518_all{nAgonist, 1}, 0, 2, 'omitnan')./sqrt(sum(~isnan(CR518_all{nAgonist, 1}), 2))         ;
    CR518_918_SEM_temp  = std(CR518_918_all{nAgonist, 1}, 0, 2, 'omitnan')./sqrt(sum(~isnan(CR518_918_all{nAgonist, 1}), 2)) ;
    GCaMP_SEM_temp      = std(GCaMP_all{nAgonist, 1}, 0, 2, 'omitnan')./sqrt(sum(~isnan(GCaMP_all{nAgonist, 1}), 2))         ;
    CR918_SEM_temp      = std(CR918_all{nAgonist, 1}, 0, 2, 'omitnan')./sqrt(sum(~isnan(CR918_all{nAgonist, 1}), 2))         ;

    % CR518
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

    % CR518-918
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

    % CR918
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

    % Find axes limits across receptors for each ligand/plot:
    maxDose_allReceptors = max([doseTempCR518(maxDose_CR518(nAgonist))
        doseTemp_518_918(maxDose_CR518_918(nAgonist))
        doseTempCR918(maxDose_CR918(nAgonist))]) ;

    minDose_allReceptors = min([doseTempCR518(1)
        doseTemp_518_918(1)
        doseTempCR918(1)]) ;
    
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

lgd = legend([h1 h2 h3], 'Orientation', 'horizontal', 'FontSize', 14);
lgd.Layout.Tile = 'south';  % or 'north', 'east', 'west'

set(gcf, 'Position',  [500, 500, 1500, 800]) % % [pos, pos, width, height]




%% MAX AMPLITUD HEATMAP
receptors = {'CR518', 'CR518-918', 'CR918'} ;

% Use maxDose for each agonist to extract values at max concentration
stackMeans = zeros(length(ligandList), 3) ;
stackSEMs  = zeros(length(ligandList), 3) ;

for nligand = 1:length(ligandList)
    % idx = maxDose_CR518(nligand);

    stackMeans(nligand,1) = mean(CR518_all{nligand}(maxDose_CR518(nligand),:), 'omitnan');
    stackMeans(nligand,2) = mean(CR518_918_all{nligand}(maxDose_CR518_918(nligand),:), 'omitnan');
    stackMeans(nligand,3) = mean(GCaMP_all{nligand}(maxDose_CR918(nligand),:), 'omitnan');

    stackSEMs(nligand,1) = std(CR518_all{nligand}(maxDose_CR518(nligand),:), 0, 2, 'omitnan') ./ sqrt(sum(~isnan(CR518_all{nligand}(maxDose_CR518(nligand),:))))   ;
    stackSEMs(nligand,2) = std(CR518_918_all{nligand}(maxDose_CR518_918(nligand),:), 0, 2, 'omitnan') ./ sqrt(sum(~isnan(CR518_918_all{nligand}(maxDose_CR518_918(nligand),:))))    ;
    stackSEMs(nligand,3) = std(GCaMP_all{nligand}(maxDose_CR918(nligand),:), 0, 2, 'omitnan') ./ sqrt(sum(~isnan(GCaMP_all{nligand}(maxDose_CR918(nligand),:))))             ;
end

figure
h = heatmap(ligandList, receptors, stackMeans', ...
    'Colormap', parula, ...
    'MissingDataColor', 'w', ...
    'MissingDataLabel', 'N/A');

title('Peak response (vehicle folds)')
h.CellLabelColor = 'none';  % Hide numbers
set(gca, 'FontSize', 12, 'FontName', 'Arial');
set(gcf, 'Position',  [500, 500, 500, 200]) % % [pos, pos, width, height]


%% EC50 HEATMAP
% Updated EC50 values (CR518 in Column1, CR518_918 in Column2)
% ENTER VALUES MANUALLY

EC50_valuesAll = nan(length(ligandList), 3); % preallocate with NaN
for nVal = 1:length(ligandList)
    % use conditional assignment with isempty
    if ~isempty(EC50(nVal).CR518)
        EC50_valuesAll(nVal,1) = EC50(nVal).CR518;
    end
    if ~isempty(EC50(nVal).CR518_918)
        EC50_valuesAll(nVal,2) = EC50(nVal).CR518_918;
    end
    if ~isempty(EC50(nVal).CR918)
        EC50_valuesAll(nVal,3) = EC50(nVal).CR918;
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
set(gcf, 'Position',  [500, 500, 500, 200]) % % [pos, pos, width, height]





