function [mutant_counts] = load_mutant_counts(sorted_reads_filename, path_to_save_output)
% loads the counts of all mutants, using the mutant identity as a binary
% vector to determine the position within the vector (all 0's is the first
% mutant, 000...01 is the second mutant, 000...10 third, so on, and all 1's is the last mutant)
   
   [unique_mutants, unique_counts]=calculate_mutant_counts(sorted_reads_filename, path_to_save_output); 
    
    mutant_index = 1+bi2de(unique_mutants); % add 1 so that all 0's is the first mutant
    mutant_counts = zeros(1024,1);
    mutant_counts(mutant_index) = unique_counts;
end