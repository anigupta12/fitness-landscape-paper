function [unique_mutants, unique_counts] = calculate_mutant_counts(sorted_reads_filename, path_to_save_output)

    % Load the data
    [seq1, seq2, counts] = load_sorted_sequences(sorted_reads_filename);

    % Load all the engineered positions in the library
    [seq1_lib_pos, seq1_lib_WT, seq1_lib_mut, seq2_lib_pos, seq2_lib_WT, seq2_lib_mut] = define_mutation_positions();
    [seq1_neut1_pos, seq1_neut1_WT, seq1_neut1_mut, seq1_neut2_pos, seq1_neut2_WT, seq1_neut2_mut, seq2_neut1_pos, seq2_neut1_WT, seq2_neut1_mut, seq2_neut2_pos, seq2_neut2_WT, seq2_neut2_mut] = define_neutral_mutation_positions();
    [seq1_rest_pos, ~, seq1_full_WT, seq2_rest_pos, ~, seq2_full_WT] = define_WT_nonEngineered_positions();

    %% Removing all the reads which have mismatch > 1 with at the non-engineered positions

    [mismatch_seq1, mismatch_per_read_all_seq1, mismatch_per_read_seq1, total_mismatch_seq1] = reads_mismatch_with_WT(seq1, seq1_full_WT, seq1_rest_pos);
    [mismatch_seq2, mismatch_per_read_all_seq2, mismatch_per_read_seq2, total_mismatch_seq2] = reads_mismatch_with_WT(seq2, seq2_full_WT, seq2_rest_pos);

    [seq1_filtered, seq2_filtered, counts_filtered, total_reads_original, total_reads_filtered, mismatch_per_read_union, mismatch_per_read_seq1, mismatch_per_read_seq2] = remove_reads_mismatch_with_WT(seq1, seq2, counts, mismatch_seq1, mismatch_seq2);

% If you want to turn off the filter then uncomment the 3 lines below
%  seq1_filtered = seq1;
%  seq2_filtered = seq2;
%  counts_filtered= counts;
    %% Parsing mutations 

    % Creating a mutant_identity matrix with 0's if there is no mutation at the
    % given position and 1 otherwise
    N = length(counts_filtered);
    for i = 1:N
        % CHANGE 'parse_neutral_mutations' below to
        % 'parse_neutral_mutations_for_neutral_effects_first5sites' for
        % effects of neutral mutations
        [mutant_identity(i,:), neut_mutant_identity(i,:)] = parse_neutral_mutations( seq1_filtered(i,:), seq2_filtered(i,:) );
    end

    %% Removing sequences that don't match WT or engineered mutations

    % Taking into account only the library mutation position (have vals = 2)
    bad_mutants = (mutant_identity == 2);
    bad_mutants_sum = sum(bad_mutants, 2);

    percent_bad = sum(bad_mutants_sum~=0)/length(bad_mutants_sum);

    %%%%----%%%%%
    % Taking into account the neutral mutation position also (have vals = -1)
    % Removing all reads that don't match neutral as well 

    bad_mutants_neutral = (mutant_identity == -1 );
    bad_mutants_total = bad_mutants + bad_mutants_neutral;
    bad_mutants_total_sum = sum(bad_mutants_total, 2);

    % Filtering the counts vector according to the filtered ("good") mutants
    good_mutants_indices = find(bad_mutants_total_sum == 0);
    good_counts = counts_filtered(good_mutants_indices,:); 
    good_mutants = mutant_identity(good_mutants_indices,:);
    
    % 
    neut_mutant_identity = mutant_identity(good_mutants_indices,:);
    
    percent_bad_neutral = sum(bad_mutants_total_sum~=0)/length(bad_mutants_total_sum);

    [bad_mutants_neutral_row, bad_mutants_neutral_column] = find(bad_mutants_total);
    bad_mutants_neutral_indices = [bad_mutants_neutral_row bad_mutants_neutral_column];
    bad_mutants_neutral_indices_sorted = sortrows(bad_mutants_neutral_indices); 


    %% Make the library unique, adding up redundant counts

    [unique_mutants, ui, mi] = unique(good_mutants, 'rows');
    unique_counts = zeros(length(ui),1);

    for i = 1:length(ui) % number of unique mutants
        redundant_indices = (mi == i);
        unique_counts(i) = sum(good_counts(redundant_indices));
    end

    %% NOT REQUIRED or CONNECTED to OTHER PARTS of CODE

    % To sort mutant matrixc and counts according to decreasing order of counts
    sort_counts = true;
    if sort_counts
        [unique_counts_sort, order] = sort(unique_counts,'descend');
        unique_mutants_sort = unique_mutants(order,:);
    end

    %% Saving entire workspace and output of each files

    formatSpec = '%s%s%s.mat';
    save_workspace_name = sprintf(formatSpec, path_to_save_output ,'workspace_parsed_mutant_', strrep(sorted_reads_filename, '/', '-')); 
    save(save_workspace_name);
    save_unique_mutant_name = sprintf(formatSpec, path_to_save_output ,'output_parsed_mutant_', strrep(sorted_reads_filename, '/', '-')); 
    save(save_unique_mutant_name,  'unique_mutants', 'unique_counts');

    %%
    % IF THERE IS ALREADY AN EXISTING OUTPUT, THEN LOAD THAT
    %     else
    %         % load the data from files 
    %         load(saved_path_parse_mutation,  'unique_mutants', 'unique_counts');
    %     end


end