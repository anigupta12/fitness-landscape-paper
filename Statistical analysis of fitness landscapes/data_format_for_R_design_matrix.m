%% Load baseData_with_selectionRate_new.mat
% Input Format -> treatment1 = [ genotypeID fitness]

treatment1 = [selectionRates_WT{1} selectionRates_WT{3}];
treatment2 = [selectionRates_malT{1} selectionRates_malT{3}];

[commonID, commonID_fitness1, commonID_fitness2] = compare_fitness_common_mutants(treatment1, treatment2);

%%
[~, malTgenotype_in_606_ID] = ismember(commonID, selectionRates_WT{1});
[~, WT_genotype_in_malT_ID] = ismember(commonID, selectionRates_malT{1});

% common_genotype_binary_ID = selectionRates_WT{2};
% common_genotype_binary_ID = common_genotype_binary_ID(malTgenotype_in_606_ID);
% 
% common_genotype_ID = selectionRates_WT{1};
% common_genotype_ID = common_genotype_ID(malTgenotype_in_606_ID);

common_genotype_fitness_WT = selectionRate_wGenotypeID_WT(malTgenotype_in_606_ID,:);

common_genotype_fitness_malT = selectionRate_wGenotypeID_malT(WT_genotype_in_malT_ID,:);

%%
design_matrix_fitness_wt = [];
design_matrix_genotypeID_wt = [];
design_matrix_genotype_wt = [];
design_matrix_host_wt = [];
design_matrix_total_mutations_wt = [];
for i = 1:length(common_genotype_fitness_WT)
    for j = 2:5
        if isnan(common_genotype_fitness_WT(i,j)) == 0
            design_matrix_fitness_wt = [design_matrix_fitness_wt; common_genotype_fitness_WT(i,j)];
            design_matrix_genotypeID_wt = [design_matrix_genotypeID_wt; common_genotype_fitness_WT(i,1)];
            design_matrix_genotype_wt = [design_matrix_genotype_wt; de2bi(common_genotype_fitness_WT(i,1),10)];
            design_matrix_total_mutations_wt = [design_matrix_total_mutations_wt; sum(de2bi(common_genotype_fitness_WT(i,1),10))];
            design_matrix_host_wt = [design_matrix_host_wt; 'wt'];
        end
    end
end
%%
design_matrix_fitness_malT = [];
design_matrix_genotypeID_malT = [];
design_matrix_genotype_malT = [];
design_matrix_host_malT = [];
design_matrix_total_mutations_malT = [];
for i = 1:length(common_genotype_fitness_malT)
    for j = 2:5
        if isnan(common_genotype_fitness_malT(i,j)) == 0
            design_matrix_fitness_malT = [design_matrix_fitness_malT; common_genotype_fitness_malT(i,j)];
            design_matrix_genotypeID_malT = [design_matrix_genotypeID_malT; common_genotype_fitness_malT(i,1)];   
            design_matrix_genotype_malT = [design_matrix_genotype_malT; de2bi(common_genotype_fitness_malT(i,1),10)];
            design_matrix_total_mutations_malT = [design_matrix_total_mutations_malT; sum(de2bi(common_genotype_fitness_malT(i,1),10))];
            design_matrix_host_malT = [design_matrix_host_malT; 'malT'];
        end
    end
end
%%
save('design_matrix_for_R_FINAL','design_matrix_fitness_wt','design_matrix_genotypeID_wt','design_matrix_genotype_wt','design_matrix_host_wt','design_matrix_fitness_malT','design_matrix_genotypeID_malT','design_matrix_genotype_malT','design_matrix_host_malT')
%%
design_matrix_fitness_mean_wt = commonID_fitness1;
design_matrix_fitness_mean_malT = commonID_fitness2;
save('design_matrix_for_R_FINAL_v2','design_matrix_fitness_mean_wt','design_matrix_total_mutations_wt','design_matrix_fitness_mean_malT','design_matrix_total_mutations_malT')
%%
csvwrite('606.csv',design_matrix_genotype_wt);
csvwrite('malT.csv',design_matrix_genotype_malT)
