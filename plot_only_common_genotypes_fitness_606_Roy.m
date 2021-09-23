treatment1 = [selectionRates_WT{1} selectionRates_WT{3}];
treatment2 = [selectionRates_malT{1} selectionRates_malT{3}];
[commonID, commonID_fitness1, ~] = compare_fitness_common_mutants(treatment1, treatment2);

%% 606
[~, index] = ismember(commonID, selectionRates_WT{1});
selectionRates_606_common = {selectionRates_WT{1}(index), selectionRates_WT{2}(index), selectionRates_WT{3}(index)};
%%
%selectionRate_wGenotypeID = selectionRate_wGenotypeID_WT;
[localpeak, neighbourMatrix, neighbourCount, localPeakCount, peakData, localPeakID, graphData] = networkGraph(selectionRates_606_common{1}, selectionRates_606_common{2}, selectionRates_606_common{3});
%%
l = length(selectionRates_606_common{:,1});

% Plotting Netowork Graph
[u] = networkGraphPlot(graphData, l, 'REL606');

%% ecc4

[~, index] = ismember(commonID, selectionRates_malT{1});
selectionRates_ecc4_common = {selectionRates_malT{1}(index), selectionRates_malT{2}(index), selectionRates_malT{3}(index)};
%%
%selectionRate_wGenotypeID = selectionRate_wGenotypeID_WT;
[localpeak, neighbourMatrix, neighbourCount, localPeakCount, peakData, localPeakID, graphData] = networkGraph(selectionRates_ecc4_common{1}, selectionRates_ecc4_common{2}, selectionRates_ecc4_common{3});
 %%
l = length(selectionRates_ecc4_common{:,1});

% Plotting Netowork Graph
[u] = networkGraphPlot(graphData, l, 'EcC4');
