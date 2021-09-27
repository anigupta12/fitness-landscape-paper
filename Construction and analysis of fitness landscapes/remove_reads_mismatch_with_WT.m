%% To remove sequencing reads which do not match with WT at non-engineered position for more than 1 bp

% Keeping all the reads which differed from WT sequence by ONLY 1 bp

function [seq1_filtered, seq2_filtered, counts_filtered, total_reads_original, total_reads_filtered, mismatch_per_read_union, mismatch_per_read_seq1, mismatch_per_read_seq2] = remove_reads_mismatch_with_WT(seq1, seq2, counts, mismatch_seq1, mismatch_seq2)

mismatch_per_read_seq1 = sum(mismatch_seq1,2); % total no. of reads x 1
mismatch_per_read_seq2 = sum(mismatch_seq2,2);

seq1_filtered = seq1;
seq2_filtered = seq2;
counts_filtered = counts;

% Taking union of reads where there is mismatch in either seq1 or seq2, as
% we need to throw away that whole read
mismatch_per_read_union = zeros(length(mismatch_per_read_seq1),1);
mismatch_per_read_union(any(mismatch_per_read_seq1>1,2),:)=1;
mismatch_per_read_union(any(mismatch_per_read_seq2>1,2),:)=1;

seq1_filtered(any(mismatch_per_read_union>0,2),:)=[];
seq2_filtered(any(mismatch_per_read_union>0,2),:)=[];
counts_filtered(any(mismatch_per_read_union>0,2),:)=[];

total_reads_original = length(seq1);  % could be seq2 or seq1
total_reads_filtered = total_reads_original - length(seq1_filtered);


end