%% Used to calculate and plot number of genotypes per n number of mutations in the final fitnness landscape
% Called in plot_distribution_of_initial_counts
function [genotype_count_per_mutation] = distribution_genotypes(genotype, L)

% Inputs are - 
% 1) a cell array like '100000000' ; '1111000010' etc.
% 2) L -> no. of loci like 10 in this case

% To convert from cell array to a matrix for bitxor operation


N = length(genotype);
ancestor_uint8 = repmat(uint8('0000000000'),N,1);
binary_genotypeID_uint8_temp = uint8(genotype{1});


for i=2:length(genotype) 
    temp = uint8(genotype{i});
    binary_genotypeID_uint8_temp = [binary_genotypeID_uint8_temp; temp];
end

binary_genotypeID_uint8 = bitxor(binary_genotypeID_uint8_temp,ancestor_uint8);

% To arrange all the genotypes in increasing number of mutations
number_mutation=zeros(N,1);
ancestor=zeros(1,L);
for i=1:N
    ancestor = binary_genotypeID_uint8(1,:);
    number_mutation(i,1) = sum(bitxor(ancestor,binary_genotypeID_uint8(i,:)));
end

[number_mutation_sorted,b] = sort(number_mutation);
binary_genotypeID_uint8_sorted = binary_genotypeID_uint8(b,:);
%neighbour_count_matrix_sorted = neighbour_count_matrix(b,:);

% To get the boundaries for change in number of mutation
boundary = zeros(L,1);

genotype_count_per_mutation = zeros(L+1,1);
genotype_count_per_mutation(1,1) = 1;
genotype_count_per_mutation(L+1,1) = 1;

for i=1:L
   [temp, ~] = find(number_mutation_sorted==i);
   boundary(i,1) = temp(1);
   genotype_count_per_mutation(i+1,1) = length(temp);
end

end