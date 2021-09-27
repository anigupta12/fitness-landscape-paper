% already have selection rates as a 1x3 cell for wt and malT

%% For WT (or 606)
selectionRates = selectionRates_WT;
selectionRate_wGenotypeID = selectionRate_wGenotypeID_WT;
[localpeak, neighbourMatrix, neighbourCount, localPeakCount, peakData, localPeakID, graphData] = networkGraph(selectionRates{1}, selectionRates{2}, selectionRates{3});

% Creating a structure of peak data with selection coefficient (fitness) and neighbours for input to ANOVA function
[peakStruct] = localPeakAnalysis(10, selectionRates{1}, selectionRates{2}, selectionRate_wGenotypeID(:,2:end), localpeak);

% ANOVA with Bonferroni Correction or Tukey's Test
[anovaOutput,  validAnovaPeaks, multiCompareOutput, validPeaks, finalpassedGenotype] = anova(peakStruct, localPeakCount, 'bonferroni', 0.05);

l = length(selectionRates{:,1});
%%
% Plotting Netowork Graph
[u] = networkGraphPlot(graphData, l, 'REL606');

%% For malT (or EcC4)
selectionRates = selectionRates_malT;
selectionRate_wGenotypeID = selectionRate_wGenotypeID_malT;
[localpeak, neighbourMatrix, neighbourCount, localPeakCount, peakData, localPeakID, graphData] = networkGraph(selectionRates{1}, selectionRates{2}, selectionRates{3});

% Creating a structure of peak data with selection coefficient (fitness) and neighbours for input to ANOVA function
[peakStruct] = localPeakAnalysis(10, selectionRates{1}, selectionRates{2}, selectionRate_wGenotypeID(:,2:end), localpeak);

% ANOVA with Bonferroni Correction or Tukey's Test
[anovaOutput,  validAnovaPeaks, multiCompareOutput, validPeaks, finalpassedGenotype] = anova(peakStruct, localPeakCount, 'bonferroni', 0.05);

l = length(selectionRates{:,1});
%%
% Plotting Netowork Graph
[u] = networkGraphPlot(graphData, l, 'EcC4');

%%