%% For fitness effects of neutral mutations
site = 6;
data_neutral1 = mutant_counts_wNeutral_wFilter_final_neutral_site6a;
data_neutral2 = mutant_counts_wNeutral_wFilter_final_neutral_site6b;
%% Create count matrix for different replicates and treatment from combined data
[phageAll_neutral1, phageInitial_neutral1] = parse_mutant_counts(data_neutral1);
[phageAll_neutral2, phageInitial_neutral2] = parse_mutant_counts(data_neutral2);

%% Calculate selection rates
[ l_WT_neutral1, selectionRates_WT_neutral1, selectionRate_wGenotypeID_WT_neutral1, selectionRate_mean_wGenotypeID_WT_neutral1] = selectionRate_new(mutant_index,phageAll_neutral1, phageInitial_neutral1, 1, 4, 1);
[ l_malT_neutral1, selectionRates_malT_neutral1, selectionRate_wGenotypeID_malT_neutral1, selectionRate_mean_wGenotypeID_malT_neutral1] = selectionRate_new(mutant_index,phageAll_neutral1, phageInitial_neutral1, 5, 8, 1);

[ l_WT_neutral2, selectionRates_WT_neutral2, selectionRate_wGenotypeID_WT_neutral2, selectionRate_mean_wGenotypeID_WT_neutral2] = selectionRate_new(mutant_index,phageAll_neutral2, phageInitial_neutral2, 1, 4, 1);
[ l_malT_neutral2, selectionRates_malT_neutral2, selectionRate_wGenotypeID_malT_neutral2, selectionRate_mean_wGenotypeID_malT_neutral2] = selectionRate_new(mutant_index,phageAll_neutral2, phageInitial_neutral2, 5, 8, 1);

%% Extract selection rate of all genotypes for site i which has the mutation on (or off, change the function for recalculate selectionRates using genotype_index)

%Input site and selection rates according to the host
% selectionRates_neutral1_input = selectionRates_WT_neutral1;
% selectionRates_neutral2_input = selectionRates_WT_neutral2;
% hostType = sprintf('WT')
% 
selectionRates_neutral1_input = selectionRates_malT_neutral1;
selectionRates_neutral2_input = selectionRates_malT_neutral2;
hostType = sprintf('malT')


[genotype_index_neutral1, selectionRates_filtered_neutral1] = get_selectionRate_wMutation_site_i(selectionRates_neutral1_input, site);
[genotype_index_neutral2, selectionRates_filtered_neutral2] = get_selectionRate_wMutation_site_i(selectionRates_neutral2_input, site);

%% Getting common genotypes present in both neutral site 1 parsing and neutral site 2, and then plotting the figures
figure()
[commonID, commonID_fitness_neutral1, commonID_fitness_neutral2 ] = compare_fitness_common_mutants( [selectionRates_filtered_neutral1{1} selectionRates_filtered_neutral1{3}], [selectionRates_filtered_neutral2{1}, selectionRates_filtered_neutral2{3}] );
[p_WT, rsq_adj_WT] = scatterPlot_wLinearFit(commonID_fitness_neutral1,commonID_fitness_neutral2);
hold on
ax = gca;
%pointsize = 20;
scatter(commonID_fitness_neutral1, commonID_fitness_neutral2,'filled')%'MarkerEdgeColor', 'black');
currentLimits = [ax.XLim, ax.YLim];
newLimits = [min(currentLimits) max(currentLimits)];
ax.XLim = newLimits;
ax.YLim = newLimits;
plot(newLimits,newLimits)
title(hostType);

% To mark the genotype with all mutations off except at site i
total_sites = 10; % Total number of sites
single_mut_genotype = cellstr(num2str(zeros(1,total_sites))); % Creating '0000000000';
single_mut_genotype = regexprep(single_mut_genotype, '\W', ''); % removing spaces between 0s
single_mut_genotype{1}(site) = '1'

[present, index_present] = ismember(single_mut_genotype, selectionRates_neutral1_input{2});
if present == true
    single_mut_ID = selectionRates_neutral1_input{1}(index_present)
    [present_in_both_neutral, single_mut_index_commonID] = ismember(single_mut_ID,commonID)
    if present_in_both_neutral == true
        scatter(commonID_fitness_neutral1(single_mut_index_commonID), commonID_fitness_neutral2(single_mut_index_commonID),'r','filled');
    end
else
    disp(present)
end

%% Extracting common genotypes in 0s and 1s corresponding the commonID

common_genotype = mutant_index(commonID+1,:);
common_genotype = cellstr(num2str(common_genotype)); 
common_genotype = regexprep(common_genotype, '\W', '');

%% TO CHECK SOMETHING.. no need to run for this neutral mutation analysis
data = mutant_counts_wNeutral_wFilter_final;
[phageAll, phageInitial] = analysis_output_of_final(data);

%% Calculate selection rates
[ l_WT, selectionRates_WT, selectionRate_wGenotypeID_WT, selectionRate_mean_wGenotypeID_WT] = selectionRate_new(mutant_index,phageAll, phageInitial, 1, 4, 1);
[ l_malT, selectionRates_malT, selectionRate_wGenotypeID_malT, selectionRate_mean_wGenotypeID_malT] = selectionRate_new(mutant_index,phageAll, phageInitial, 5, 8, 1);
