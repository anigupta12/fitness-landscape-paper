function [phageAll, phageInitial] = parse_mutant_counts(data)

% Dealing with entries which didn't show up after sequencing

NaNindices=data(:)==0.1;  % Removing all the entries which are 0.1 and represent 0 actually
dataZero=data;
dataZero(NaNindices)=0;

%% Trying to create the matrix with all the mutations (checkered image) and other 4 and 4 columns
% phageInitial_a=mean(dataZero(:,29:30),2); 
% phageInitial_b=mean(dataZero(:,31:32),2); 

%phageInitial_NaN=mean(dataNaN(:,17:18),2);    % mean of 1 and NaN will be NaN

phageInitial=mean(dataZero(:,17:18),2);


phageWTa=mean(dataZero(:,1:2),2);
phageWTb=mean(dataZero(:,3:4),2);
phageWTc=mean(dataZero(:,5:6),2);
phageWTd=mean(dataZero(:,7:8),2);
phageECCa=mean(dataZero(:,9:10),2);
phageECCb=mean(dataZero(:,11:12),2);
phageECCc=mean(dataZero(:,13:14),2);
phageECCd=mean(dataZero(:,15:16),2);

phageAll = [phageWTa phageWTb phageWTc phageWTd phageECCa phageECCb phageECCc phageECCd];
%phageInitial = [phageInitial_a phageInitial_b];

% sth wrong with this
phageAll_sum=repmat(sum(phageAll),size(phageAll,1),1);  % sum of each column of phageAll and repeating that no. for n rows to make same size matrix as phageAll
phageAll_freq=phageAll./phageAll_sum; 

end
