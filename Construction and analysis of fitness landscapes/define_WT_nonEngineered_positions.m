function [seq1_rest_pos, seq1_rest_WT, seq1_full_WT, seq2_rest_pos, seq2_rest_WT, seq2_full_WT] = define_WT_nonEngineered_positions()
    
    [seq1_lib_pos, ~, ~, seq2_lib_pos, ~, ~] = define_mutation_positions();
    [seq1_neut1_pos, ~, ~, seq1_neut2_pos, ~, ~, seq2_neut1_pos, ~, ~, seq2_neut2_pos, ~, ~] = define_neutral_mutation_positions();
    
    seq1_rest_pos = 1:72;
    seq1_rest_pos = setdiff(seq1_rest_pos, seq1_lib_pos);
    seq1_rest_pos = setdiff(seq1_rest_pos, seq1_neut1_pos);
    seq1_rest_pos = setdiff(seq1_rest_pos, seq1_neut2_pos);
    
    % full WT sequence1 
    seq1_full_WT = 'CGGAAAAAATCGTCGGGGACATTGTAAAGGCGGCGAGCGCGGCTTTTCCGCGCCAGCGTGAAAGCAGTGTGG';
    
    seq1_rest_WT = seq1_full_WT;
    seq1_rest_WT(seq1_lib_pos) = 'N';
    seq1_rest_WT(seq1_neut1_pos) = 'N';
    seq1_rest_WT(seq1_neut2_pos) = 'N';
    [~, b] = find(seq1_rest_WT=='N');
    
    % WT sequence1 without engineered positions
    seq1_rest_WT(b) = [];
    
    seq2_rest_pos = 1:72;
    seq2_rest_pos = setdiff(seq2_rest_pos, seq2_lib_pos);
    seq2_rest_pos = setdiff(seq2_rest_pos, seq2_neut1_pos);
    seq2_rest_pos = setdiff(seq2_rest_pos, seq2_neut2_pos);
    
    % full WT sequence2
    seq2_full_WT = 'GCATTCGGCAGATATTCCGCCGTATACGTTTGCCAGCGATGTGCAGGTTATGGTGATTAAGAAACAGGCGCT';
    
    seq2_rest_WT = seq2_full_WT;
    seq2_rest_WT(seq2_lib_pos) = 'N';
    seq2_rest_WT(seq2_neut1_pos) = 'N';
    seq2_rest_WT(seq2_neut2_pos) = 'N';
    [~, b] = find(seq2_rest_WT=='N');
    
    % WT sequence2 without engineered positions
    seq2_rest_WT(b) = [];
    
end