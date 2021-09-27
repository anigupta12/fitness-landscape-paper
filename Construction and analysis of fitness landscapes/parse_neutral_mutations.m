function [mutant_identity, neut_mutant_identity] = parse_neutral_mutations(seq1, seq2)

[seq1_lib_pos, seq1_lib_WT, seq1_lib_mut, seq2_lib_pos, seq2_lib_WT, seq2_lib_mut] = define_mutation_positions();
[seq1_neut1_pos, seq1_neut1_WT, seq1_neut1_mut, seq1_neut2_pos, seq1_neut2_WT, seq1_neut2_mut, seq2_neut1_pos, seq2_neut1_WT, seq2_neut1_mut, seq2_neut2_pos, seq2_neut2_WT, seq2_neut2_mut] = define_neutral_mutation_positions();


size_lib_seq1 = length(seq1_lib_pos);
size_lib_seq2 = length(seq2_lib_pos);

mutant_identity = -1*ones(1,size_lib_seq1 + size_lib_seq2);
neut_mutant_identity = -1*ones(1,size_lib_seq1 + size_lib_seq2);

for i = 1:size_lib_seq1
    pos = seq1_lib_pos(i);  % eg. pos = 21
    nt = seq1(pos);           % eg. nt = s1(21), say G
    pos_neut1 = seq1_neut1_pos(i);
    nt_neut1 = seq1(pos_neut1);           
    pos_neut2 = seq1_neut2_pos(i);
    nt_neut2 = seq1(pos_neut2);           

    if nt == seq1_lib_WT(i)
        if nt_neut1 == seq1_neut1_WT(i) || nt_neut2 == seq1_neut2_WT(i)
        %if nt_neut2 == seq1_neut2_WT(i)
            mutant_identity(i) = 0;
            neut_mutant_identity(i) = 0;
        end
    elseif nt == seq1_lib_mut(i)
        if nt_neut1 == seq1_neut1_mut(i) || nt_neut2 == seq1_neut2_mut(i)
        %if nt_neut2 == seq1_neut2_mut(i)
            mutant_identity(i) = 1;
            neut_mutant_identity(i) = 1;
        end

    else
        mutant_identity(i) = 2;
    end


end

for i = 1:size_lib_seq2
    pos = seq2_lib_pos(i);
    nt = seq2(pos);
    pos_neut1 = seq2_neut1_pos(i);
    nt_neut1 = seq2(pos_neut1);           
    pos_neut2 = seq2_neut2_pos(i);
    nt_neut2 = seq2(pos_neut2);           

    if nt == seq2_lib_WT(i)  
        if nt_neut1 == seq2_neut1_WT(i) || nt_neut2 == seq2_neut2_WT(i)
        %if nt_neut2 == seq2_neut2_WT(i)
            mutant_identity(size_lib_seq1+i) = 0;
            neut_mutant_identity(size_lib_seq1+i) = 0;
        end
    elseif nt == seq2_lib_mut(i)
        if nt_neut1 == seq2_neut1_mut(i) || nt_neut2 == seq2_neut2_mut(i)
        %if nt_neut2 == seq2_neut2_mut(i)
            mutant_identity(size_lib_seq1+i) = 1;
            neut_mutant_identity(size_lib_seq1+i) = 1;
        end
    else
        mutant_identity(size_lib_seq1+i) = 2;
    end
end

end
