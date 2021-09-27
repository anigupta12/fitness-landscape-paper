function [ l, selectionRates, selectionRate_wGenotypeID, selectionRate_mean_wGenotypeID] = selectionRate_new(mutant_index, phageAll_Final, phageInitial_NaN, a, b, wt)

% Output- 
% a. l (scalar)- total no. of valid detected genotypes
% b. selectionRates (1x3 cell}- 1st column is genotypeID (decimal), 2nd- genotype (binary), 3rd- mean selection rate of genotype
% c. selectionRate_wGenotyeID (lx(1+r) matrix)- selection rates of all the detected genotypes; 1st column is genotypeID in decimal; r is the no. of repliactes for the given treatment
% d. selectionRate_mean_wGenotypeID (l x (1+1) matrix)- mean of selection rates of genptypes over all the replicates; 1st column is genotypeID in decimal


% Input-
% a. phageAll (2^n x m matrix)- this contains counts for different replicates for all the treatments; n is total no. of point mutations and m = total no. of replicates x no. of treatments
% b. phageInitial_NaN (2^n x 1 matrix)- Initial counts of all the genotypes before starting the experiment with Nan for undetected genotypes
% c. a (scalar)- starting column containing counts for the given treatment; b-a+1 = total no. of replicates
% d. b (scalar)- last column containing counts of the given treatment
% e. wt (scalar)- index/row containging Initial counts of wt phage, i.e. reference for calculating selection rate

phageInitialWT = phageInitial_NaN(wt);
phageFinalWT = phageAll_Final(1,a:b); %matrix of 1x4 for 4 replicates
logWT = log(phageFinalWT./phageInitialWT); %matrix of 1024x4
phageFinal = phageAll_Final(:,a:b); %matrix of 1024x4 for 4 replicates

N = length(phageAll_Final); % total no. of genotypes in phageAll matrix, eg. 1024 if there are 10 sites of mutation
phageInitial_NaN = repmat(phageInitial_NaN,1,b-a+1); %creating matrix of 1024x4 from 1024x1 because initial seed was same in all the replicates
logWT = repmat(logWT,N,1); %creating matrix of 2^10=1024 x 4 from 1x4 because ln(finalWT/initialWT) will be same for all genotypes for different replicates

selectionRate = log(phageFinal./phageInitial_NaN)-logWT; %matrix of 1024x4 for selection rates

% Replacing all the inf values with NaN
selectionRate_NaN = selectionRate;
selectionRate_NaN(selectionRate == inf) = nan;
selectionRate_NaN(selectionRate_NaN == -inf) = nan;

% selectionRate with genotypeID (in decimal)
selectionRate_wGenotypeID_NaN = [transp(0:N-1)  selectionRate_NaN];
% selection rate with genotypeID with only detected genotypes  
selectionRate_wGenotypeID = selectionRate_wGenotypeID_NaN(~all(isnan(selectionRate_wGenotypeID_NaN(:,2:end)),2),:);
% total no. of detected genotypes
l = length(selectionRate_wGenotypeID);
% mean of selection rates across replicates with NaN genotypes
selectionRate_mean_NaN =  nanmean(selectionRate_NaN,2);
% mean of selection rates across replicates with only detected genotypes
selectionRate_mean_wGenotypeID = [selectionRate_wGenotypeID(:,1) nanmean(selectionRate_wGenotypeID(:,2:end),2)];

genotypeID = selectionRate_wGenotypeID(:,1); % decimal representation of genotypes

%genotype = cellstr(dec2bin(genotypeID)); % binary representation of genotypes; a lx1 cell; where l is the no. of detected genotypes 
genotype = cellstr(num2str(mutant_index));
genotype = regexprep(genotype, '\W', ''); % all 1024 genotypes present
genotype = genotype(genotypeID+1); % only the genotypes present whose selection rate is available

selectionRates = {genotypeID, genotype, selectionRate_mean_wGenotypeID(:,2)};

end
