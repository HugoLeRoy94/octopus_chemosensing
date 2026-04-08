
% cd to folder containing '20251017_HiPlexResults.mat'

clear
close all
clc

load('20251017_HiPlexResults.mat')


%% RASTER PLOT WITH ALL CELLS
% In 'binaryTable_Stack' 1= presence, 0= absence 

figure, imagesc(binaryTable_Stack)
xlabel('Cell number')
set(gca,'YTick',1:size(CRnames, 1), 'YTickLabel',CRnames, 'TickLabelInterpreter','none')
set(gcf, 'Position',  [250, 250, 900, 400]) % width by height in pixels
colormap(flipud(gray)) % cool summer winter autumn
set(gca,'TickLength',[0 0])
set(gca, 'FontSize', 12)


%% HISTOGRAM WITH NUMBER OF COEXPRESSED SUBUNITS
coexpMat= sum(binaryTable_Stack, 1) ;
coexpMat(coexpMat == 0) = []        ; 

[counts, edges] = histcounts(coexpMat, 'BinMethod','integers') ;
centers = edges(1:end-1) + 0.5; % bin centers

figure, bar(centers, counts, 0.7, 'FaceColor',[0.6 0.6 0.6]) 
ylabel('Number of cells')
xlabel('Number of coexpressed receptors')
set(gca, 'FontSize', 12)


%% RELATIVE ABUNDANCE OF EACH RECEPTOR IN THE POPULATION
CR_counts = sum(binaryTable_Stack, 2) ;
figure, bar(CR_counts, 0.7, 'FaceColor',[0.6 0.6 0.6])
ylabel('Number of cells')
set(gca,'XTick',1:size(CRnames, 1), 'XTickLabel',CRnames, 'TickLabelInterpreter','none')
set(gca, 'FontSize', 12)
set(gcf, 'Position',  [250, 250, 500, 250]) % width by height in pixels
set(gca,'TickLength',[0 0])



%% PAIRS LIKELIHOODS

%  CO-OCCURRENCE PROBABILITY
pairComparison = [  ] ; 
pairComparison = binaryTable_Stack * binaryTable_Stack.' ;
pairProb = pairComparison / size(binaryTable_Stack, 2) ;

figure, imagesc(pairProb)
xlabel('Cell number')
ylabel('Cell number')
title ('Co-occurrence probability')
set(gca,'XTick',1:size(CRnames, 1), 'XTickLabel',CRnames, 'TickLabelInterpreter','none')
set(gca,'YTick',1:size(CRnames, 1), 'YTickLabel',CRnames, 'TickLabelInterpreter','none')
set(gcf, 'Position',  [250, 250, 500, 450]) % width by height in pixels
colormap(flipud(gray)) % cool summer winter autumn
set(gca,'TickLength',[0 0])
set(gca, 'FontSize', 12)
xtickangle(90)


%%  CONDITIONAL PROBABILITY
numCells_perCR = diag(pairComparison); % cells expressing each receptor
condProb = pairComparison ./ numCells_perCR;

figure, imagesc(condProb)
xlabel('Receptor A')
ylabel('Receptor B')
title ('Conditional probability (P(A|B))')
set(gca,'XTick',1:size(CRnames, 1), 'XTickLabel',CRnames, 'TickLabelInterpreter','none')
set(gca,'YTick',1:size(CRnames, 1), 'YTickLabel',CRnames, 'TickLabelInterpreter','none')
set(gcf, 'Position',  [250, 250, 500, 450]) % width by height in pixels
colormap(flipud(gray)) % cool summer winter autumn
set(gca,'TickLength',[0 0])
set(gca, 'FontSize', 12)
xtickangle(90)



%% NUMBER OF UNIQUE COMBINATIONS
M_cells = binaryTable_Stack';   % size ~2000 x 26

[uniquePatterns, ~, idxUnique] = unique(M_cells, 'row', 'stable');
uniqueReceptors = size(uniquePatterns, 1);
counts = accumarray(idxUnique, 1);

fprintf('Number of unique expression patterns: %d\n', uniqueReceptors);

% Most common patterns
[sortedCounts, order] = sort(counts, 'descend');
topPatt = min(100,uniqueReceptors);
disp('Top patterns:');
disp(sortedCounts(1:topPatt));

for nPatt = 1:topPatt
    patternIdx = order(nPatt);
    genesOn = CRnames(uniquePatterns(patternIdx,:) == 1);
    fprintf('\nPattern %d (appears %d times):\n', nPatt, sortedCounts(nPatt));
    disp(genesOn');
end

%% BATCH ANALYSIS-MOVING WINDOWS
M_cells = binaryTable_Stack';   % size ~2000 x 26
nCells = size(M_cells,1);

batchEdge = 10:10:nCells;
nBatches = length(batchEdge) ;
uniquePerBatch = zeros(1, nBatches);
startIdx = 1;

for nBatch = 1:nBatches
    stopIdx = min(startIdx + batchEdge(nBatch) - 1, nCells);
    if startIdx > nCells
        break; % last batch
    end   
    subSet_temp = M_cells(startIdx:stopIdx, :) ;  % SUBSET WINDOW TEMP
    uniquePerBatch(nBatch) = size(unique(subSet_temp, 'rows'), 1); % UNIQUE PATTERNS
    %startIdx = stopIdx + 1; % MOVES THE WINDOW
end

figure;
bar(uniquePerBatch)
xlabel('Window number')
ylabel('Unique # of patterns in batch')
%xlim([0 20])
xticks(1:nBatches)
title('Overlapping batches')
box off

% lower number of x-axis labels
numLabels = 10;  % label ticks number
idxLabels = round(linspace(1, nBatches, numLabels));  % even spaces
xticks(idxLabels)
xticklabels(arrayfun(@num2str, idxLabels, 'UniformOutput', false))


% clearvars -except binaryTable_Stack CRnames


