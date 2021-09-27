
data = mutant_counts_wNeutral_wFilter;
%% Map mutant_counts data to different treatments
[phageAll, phageInitial] = parse_mutant_counts(data);


%% Get selection rates 
[ l_WT, selectionRates_WT, selectionRate_wGenotypeID_WT, selectionRate_mean_wGenotypeID_WT] = selectionRate_new(mutant_index,phageAll, phageInitial, 1, 4, 1);
[ l_malT, selectionRates_malT, selectionRate_wGenotypeID_malT, selectionRate_mean_wGenotypeID_malT] = selectionRate_new(mutant_index,phageAll, phageInitial, 5, 8, 1);

%% Get selection rates for neutral 1 vs neutral 2 fitness effect


[ l_WT_neutral1, selectionRates_WT_neutral1, selectionRate_wGenotypeID_WT_neutral1, selectionRate_mean_wGenotypeID_WT_neutral1] = selectionRate_new(mutant_index,phageAll_neutral1, phageInitial_neutral1, 1, 4, 1);
[ l_malT_neutral1, selectionRates_malT_neutral1, selectionRate_wGenotypeID_malT_neutral1, selectionRate_mean_wGenotypeID_malT_neutral1] = selectionRate_new(mutant_index,phageAll_neutral1, phageInitial_neutral1, 5, 8, 1);

[ l_WT_neutral2, selectionRates_WT_neutral2, selectionRate_wGenotypeID_WT_neutral2, selectionRate_mean_wGenotypeID_WT_neutral2] = selectionRate_new(mutant_index,phageAll_neutral2, phageInitial_neutral2, 1, 4, 1);
[ l_malT_neutral2, selectionRates_malT_neutral2, selectionRate_wGenotypeID_malT_neutral2, selectionRate_mean_wGenotypeID_malT_neutral2] = selectionRate_new(mutant_index,phageAll_neutral2, phageInitial_neutral2, 5, 8, 1);


%%

%[~,idx] = ismember(selectionRate_mean_wGenotypeID_WT(:,1)==selectionRate_mean_wGenotypeID_WT_neutral(:,1));

%[sharedVals,idxsIntoA] = intersect(selectionRate_mean_wGenotypeID_WT(:,1),selectionRate_mean_wGenotypeID_WT_neutral(:,1),'stable')

common_genotypeID_WT = selectionRate_mean_wGenotypeID_WT(:,1);
common_genotypeID_WT  = common_genotypeID_WT (ismember(common_genotypeID_WT,selectionRate_mean_wGenotypeID_WT_neutral(:,1) ) );

Z1 = ismember(selectionRate_mean_wGenotypeID_WT(:,1),common_genotypeID_WT);
Z2 = selectionRate_mean_wGenotypeID_WT(:,2);
common_genotype_selection_rate_WT = Z2(Z1>0);

Z1 = ismember(selectionRate_mean_wGenotypeID_WT_neutral(:,1),common_genotypeID_WT);
Z2 = selectionRate_mean_wGenotypeID_WT_neutral(:,2);
common_genotype_selection_rate_WT_neutral = Z2(Z1>0);
%%
figure(1)
scatter(common_genotype_selection_rate_WT,common_genotype_selection_rate_WT_neutral)
hold on
x=[-3:11];
y=x;
plot(x,y)
title('606')
xlabel('w/o neutral mutation parsing')
ylabel('with neutral mutation parsing')

%%
common_genotypeID_malT = selectionRate_mean_wGenotypeID_malT(:,1);
common_genotypeID_malT  = common_genotypeID_malT (ismember(common_genotypeID_malT,selectionRate_mean_wGenotypeID_malT_neutral(:,1) ) );

Z1 = ismember(selectionRate_mean_wGenotypeID_malT(:,1),common_genotypeID_malT);
Z2 = selectionRate_mean_wGenotypeID_malT(:,2);
common_genotype_selection_rate_malT = Z2(Z1>0);

Z1 = ismember(selectionRate_mean_wGenotypeID_malT_neutral(:,1),common_genotypeID_malT);
Z2 = selectionRate_mean_wGenotypeID_malT_neutral(:,2);
common_genotype_selection_rate_malT_neutral = Z2(Z1>0);
%%
figure(2)
scatter(common_genotype_selection_rate_malT,common_genotype_selection_rate_malT_neutral)
hold on
x=[-3:11];
y=x;
plot(x,y)
title('ECC4')
xlabel('w/o neutral mutation parsing')
ylabel('with neutral mutation parsing')

%%
figure(3)

%barX = char(['606 '; 'ECC4']);
%barX = [1 1 ;2 2];

barY = [mean(logWT_606_wFilter_woNeutral) mean(logWT_606_wNeutral) ; mean(logWT_ECC4_wFilter_woNeutral) mean(logWT_ECC4_wNeutral)];

xTickLabels = {'606' 'ECC4'};

A = [ 1 1; 2 2]
barX = categorical(A,[1 2],{'606' 'ECC4'})

stdDev = [std(logWT_606_wFilter_woNeutral) std(logWT_606_wNeutral) ; std(logWT_ECC4_wFilter_woNeutral) std(logWT_ECC4_wNeutral)];

bar(barY);
set(gca,'XTickLabel',{'606' 'ECC4'})
%hold on
%errorbar(barY,stdDev, '.')
%%

