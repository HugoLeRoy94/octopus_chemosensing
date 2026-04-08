%% HETEROMER DATA ANALYSIS BELOW
% clearvars -except vectorAll CR_names moleculesID

% NAN AND ZEROS IN 'vectorAll' MEANS NO RESPONSE 

pearsonCorr = [ ] ; % Pearson correlation

idx_CR518 = strfind(CR_names, '000-000-CR518') ; % CR518=> BROADEST TUNING; Chnage here to test other comparisons  
CR518_homomerIDX = find(~cellfun(@isempty, idx_CR518)) ;
ref_CR518_homomer = vectorAll(CR518_homomerIDX, :)     ;

% CORRELATION CR518 AND HETEROMERS FOR ALL MOLECULES
for j = 1: size(vectorAll,1)
    vector_j = vectorAll(j, :) ;
    nonzeroIDX    = find(sum([ref_CR518_homomer; vector_j],1) ~=0) ; % nonzero elements common for both vectors
    tempCorr = corrcoef(ref_CR518_homomer(nonzeroIDX),vector_j(nonzeroIDX));
    pearsonCorr(j) = tempCorr(2, 1)  ;
end

[~, receptorCorr_IDX] = sort(pearsonCorr, 'descend')    ;

% SORT RECEPTORS BY CORRELATION TO CR518
CRnames_sorted= { } ;
for i= 1: length(pearsonCorr)
    CRnames_sorted(i) = CR_names(receptorCorr_IDX(i)) ;
end

% REMOVE NAN, NEGATIVE (DMSO) AND POSITIVE (IONOMYCIN) CONTROL CELLS
idxNaN_temp= contains(moleculesID, {'NaN','dmso', 'Ionomycin'}) ;
idxNaN= find(idxNaN_temp)           ;
vectorAll_noNaN= vectorAll          ;
vectorAll_noNaN(:, idxNaN) = [  ]   ;
moleculesID_noNaN= moleculesID      ;
moleculesID_noNaN(:, idxNaN) = [  ] ; 

vectorAll_mean= mean(vectorAll_noNaN, 1, "omitmissing") ;
[~, sortLigands_IDX]= sort(vectorAll_mean, 'descend')    ;

%% PLOTTING
% Matrix with all receptors sorted by corr to CR518 and ligands sorted by
% mean reponse
figure, imagesc(vectorAll_noNaN(receptorCorr_IDX, sortLigands_IDX)) ; % only CR518 vs all 
sortedMolecNames = { } ;
for i= 1: length(moleculesID_noNaN)
    sortedMolecNames(i) = moleculesID_noNaN(sortLigands_IDX(i)) ;
end

xticks(1:length(sortedMolecNames))  ;
yticks(1:length(CRnames_sorted))  ;
set(gca,'XTickLabel',sortedMolecNames, 'TickLabelInterpreter', 'none') ;
set(gca,'YTickLabel',CRnames_sorted, 'FontSize', 12, 'TickLabelInterpreter', 'none') ;

hold on
nx = size(vectorAll_noNaN, 2)   ;
ny = size(vectorAll_noNaN, 1)   ;
edge_x = repmat((0: nx)+ 0.5, ny+1, 1)   ;
edge_y = repmat((0: ny)+ 0.5, nx+1, 1).' ;
plot(edge_x ,edge_y, 'Color', 'k')     % vertical lines
plot(edge_x.', edge_y.','Color', 'k')  % horizontal lines
set(gcf, 'Position',  [250, 250, 1000, 800]) % width by height in pixels
colormap(jet)
set(gca,'FontSize',13), colorbar eastoutside;
ax = gca; ax.TickLength = [0 0];   % removes tick lines
title('Normalized response')
ylabel('Sorted by corr to CR518')
xlabel('Sorted by mean ligand response')

%% PLOTTING
% Matrix with all receptors with lignads sorted by mean response
figure, imagesc(vectorAll_noNaN(:, sortLigands_IDX)) ; % only CR518 vs all 

xticks(1:length(moleculesID_noNaN))  ;
yticks(1:length(CR_names))  ;
set(gca,'XTickLabel',sortedMolecNames, 'TickLabelInterpreter', 'none') ;
set(gca,'YTickLabel',CR_names, 'FontSize', 12, 'TickLabelInterpreter', 'none') ;

hold on
nx = size(vectorAll, 2)   ;
ny = size(vectorAll, 1)   ;
edge_x = repmat((0: nx)+ 0.5, ny+1, 1)   ;
edge_y = repmat((0: ny)+ 0.5, nx+1, 1).' ;
plot(edge_x ,edge_y, 'Color', 'k')     % vertical lines
plot(edge_x.', edge_y.','Color', 'k')  % horizontal lines
set(gcf, 'Position',  [250, 250, 1000, 800]) % width by height in pixels
colormap(jet)
set(gca,'FontSize',13), colorbar eastoutside;
ax = gca; ax.TickLength = [0 0];   % removes tick lines
title('Normalized response')
xlabel('Sorted by mean ligand response')


save('ligandScreen_variables.mat')


