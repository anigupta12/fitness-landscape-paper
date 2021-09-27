function [seq1_neut1_pos, seq1_neut1_WT, seq1_neut1_mut, seq1_neut2_pos, seq1_neut2_WT, seq1_neut2_mut, seq2_neut1_pos, seq2_neut1_WT, seq2_neut1_mut, seq2_neut2_pos, seq2_neut2_WT, seq2_neut2_mut] = define_neutral_mutation_positions()
    
    seq1_neut1_pos = [2 17 26 32 65];
    seq1_neut1_WT = 'GGAGC';
    seq1_neut1_mut = 'CTCAT';

    seq1_neut2_pos = [2 20 26 32 68];
    seq1_neut2_WT = 'GCAGT';
    seq1_neut2_mut = 'TTGCA';

    seq2_neut1_pos = [1 7 10 16 70];
    seq2_neut1_WT = 'GGATG';
    seq2_neut1_mut = 'ACCAA';

    seq2_neut2_pos = [1 7 10 16 70];
    seq2_neut2_WT = 'GGATG';
    seq2_neut2_mut = 'TAGCC';

end
