function [seq1_lib_pos, seq1_lib_WT, seq1_lib_mut, seq2_lib_pos, seq2_lib_WT, seq2_lib_mut] = define_mutation_positions()
    
    seq1_lib_pos = [1 21 23 31 66];
    seq1_lib_WT = 'CATCA';
    seq1_lib_mut = 'TGGTG';

    seq2_lib_pos = [2, 11, 12, 13, 72];
    seq2_lib_WT = 'CGATT';
    seq2_lib_mut = 'TAGAC';
    
end