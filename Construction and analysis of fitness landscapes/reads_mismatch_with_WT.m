 %% To compute the mismatches of the reads with the sequencing reads
 
 function [mismatch, mismatch_per_read_all, mismatch_per_read, total_mismatch] = reads_mismatch_with_WT(seq, seq_full_WT, seq_rest_pos)
    % length=no. of reads in seq 1, no. of columns = no. of non engineered (WT) positions in sequence;
    mismatch = not(seq(:,seq_rest_pos)==seq_full_WT(seq_rest_pos));   % 0 if base matches and 1 otherwise
   
   mismatch_per_site = sum(mismatch,1); % 1 x no. of non-engineered sites
 %  plot(mismatch_per_site)
   
   mismatch_per_read_all = sum(mismatch,2); % total no. of reads x 1
   
   mismatch_per_read = sum(mismatch,2);
   mismatch_per_read(any(mismatch_per_read==0,2),:)=[];  % total of no. reads which did not match completely with the WT sequence for non-engineered position only
   total_mismatch = length(mismatch_per_read);
   
     % plotting histogram
%    h = histogram(mismatch_per_read);
%    %plot(mismatch_per_read)
%    
%    figure(2)
%    % to plot the distribution of the positions vs mismatch to see if there
%     % is any bias in some base positions which often mutates
%    plot(mismatch_per_site);
   
 end