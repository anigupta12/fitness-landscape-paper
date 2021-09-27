    %% To create the final mutant_count matrix. CHECK the PARSE_NEUTRAL_MUTATION file for neutral mutation parsing
file_list = data_get_all_file_names('reads');
%file_list = file_list(2:end);

mutant_counts_wNeutral_wFilter = zeros(1024,numel(file_list));

for i=1:numel(file_list)
    sorted_reads_filename = file_list{i};
    path_to_save_output = '../Run 1';
    mutant_counts_wNeutral_wFilter(:,i) = load_mutant_counts(sorted_reads_filename, path_to_save_output); 
    i
end

mutant_counts_wNeutral_wFilter(mutant_counts_wNeutral_wFilter==0)=.1; 
mutant_index = de2bi(0:1023);

saved_path = sprintf('../Run 1');
save(saved_path)%,  'mutant_counts_neutral', 'mutant_index');

%% Some plot tests
% index = 1;
% figure(1);clf;
% plot(log10(mutant_counts_wNeutral_wFilter(:,index)));axis tight

% figure(2);clf;
% imagesc(log10(mutant_counts_wNeutral_wFilter));
% colormap(gray);

% figure(3);clf;
% plot(log10(mc(:,1)), log10(mc(:,2)), 'o')

%% %% Map mutant_counts data to different host treatments
data = mutant_counts_wNeutral_wFilter
[phageAll, phageInitial] = parse_mutant_counts(data);


%% Get selection rates 
[ l_WT, selectionRates_WT, selectionRate_wGenotypeID_WT, selectionRate_mean_wGenotypeID_WT] = selectionRate_new(mutant_index,phageAll, phageInitial, 1, 4, 1);
[ l_malT, selectionRates_malT, selectionRate_wGenotypeID_malT, selectionRate_mean_wGenotypeID_malT] = selectionRate_new(mutant_index,phageAll, phageInitial, 5, 8, 1);

%% Get selection rates for neutral 1 vs neutral 2 fitness effect


[ l_WT_neutral1, selectionRates_WT_neutral1, selectionRate_wGenotypeID_WT_neutral1, selectionRate_mean_wGenotypeID_WT_neutral1] = selectionRate_new(mutant_index,phageAll_neutral1, phageInitial_neutral1, 1, 4, 1);
[ l_malT_neutral1, selectionRates_malT_neutral1, selectionRate_wGenotypeID_malT_neutral1, selectionRate_mean_wGenotypeID_malT_neutral1] = selectionRate_new(mutant_index,phageAll_neutral1, phageInitial_neutral1, 5, 8, 1);

[ l_WT_neutral2, selectionRates_WT_neutral2, selectionRate_wGenotypeID_WT_neutral2, selectionRate_mean_wGenotypeID_WT_neutral2] = selectionRate_new(mutant_index,phageAll_neutral2, phageInitial_neutral2, 1, 4, 1);
[ l_malT_neutral2, selectionRates_malT_neutral2, selectionRate_wGenotypeID_malT_neutral2, selectionRate_mean_wGenotypeID_malT_neutral2] = selectionRate_new(mutant_index,phageAll_neutral2, phageInitial_neutral2, 5, 8, 1);

