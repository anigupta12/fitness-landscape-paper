%% To calculate and plot number of genotypes per n number of mutations in the final fitnness landscape

L = 10; % Maximum number of mutations or number of loci

genotype_count_per_mutation_WT = distribution_genotypes(selectionRates_WT{2}, L);
genotype_count_per_mutation_malT = distribution_genotypes(selectionRates_malT{2}, L);

%% Combined plots for WT and malT

fig = figure()
hold on
subplot(2,1,1)
bar1 = bar(0:L,genotype_count_per_mutation_WT );
bar1.FaceColor = [0.75 0.75 0.75];
ymax = max(genotype_count_per_mutation_WT);
ylim([0  ymax+0.1*ymax]);
text(0:length(genotype_count_per_mutation_WT)-1,genotype_count_per_mutation_WT,num2str(genotype_count_per_mutation_WT),'vert','bottom','horiz','center'); 
box off
title('ancestor', 'FontSize',14, 'FontName','Arial')
xlabel('Number of mutations', 'FontSize',12, 'FontName','Arial')
ylabel('Number of genotypes', 'FontSize',12, 'FontName','Arial')

subplot(2,1,2)
bar2 = bar(0:L,genotype_count_per_mutation_malT)
bar2.FaceColor = [0.75 0.75 0.75];
ymax = max(genotype_count_per_mutation_malT);
ylim([0  ymax+0.1*ymax]);
text(0:length(genotype_count_per_mutation_malT)-1,genotype_count_per_mutation_malT,num2str(genotype_count_per_mutation_malT),'vert','bottom','horiz','center'); 
box off
title('\it malT^{?}', 'FontSize',14, 'FontName','Arial')
xlabel('Number of mutations', 'FontSize',12, 'FontName','Arial')
ylabel('Number of genotypes', 'FontSize',12, 'FontName','Arial')

%% Plot for total number of possible genotypes with given number of mutations

total_genotypes_binary  = cellstr(dec2bin(1:2^L-1));
genotype_count_per_mutation_ideal = distribution_genotypes(total_genotypes_binary, L);

fig = figure()
bar1 = bar(0:L,genotype_count_per_mutation_ideal );
bar1.FaceColor = [0.75 0.75 0.75];
ymax = max(genotype_count_per_mutation_ideal);
ylim([0  ymax+0.1*ymax]);
text(0:length(genotype_count_per_mutation_ideal)-1,genotype_count_per_mutation_ideal,num2str(genotype_count_per_mutation_ideal),'vert','bottom','horiz','center'); 
box off
title('ancestor', 'FontSize',14, 'FontName','Arial')
xlabel('Number of mutations', 'FontSize',12, 'FontName','Arial')
ylabel('Number of genotypes', 'FontSize',12, 'FontName','Arial')


%% Independent plots for WT and malT

% figure()
% bar(0:L,genotype_count_per_mutation_WT)
% ymax = max(genotype_count_per_mutation_WT);
% ylim([0  ymax+0.1*ymax]);
% text(0:length(genotype_count_per_mutation_WT)-1,genotype_count_per_mutation_WT,num2str(genotype_count_per_mutation_WT),'vert','bottom','horiz','center'); 
% box off
% title('ancestor')
% xlabel('Number of mutations')
% ylabel('Number of genotypes')
% 
% %
% figure()
% bar(0:L,genotype_count_per_mutation_malT)
% ymax = max(genotype_count_per_mutation_malT);
% ylim([0  ymax+0.1*ymax]);
% text(0:length(genotype_count_per_mutation_malT)-1,genotype_count_per_mutation_malT,num2str(genotype_count_per_mutation_malT),'vert','bottom','horiz','center'); 
% box off
% title('malT-')
% xlabel('Number of mutations')
% ylabel('Number of genotypes')

%% CHECK CODE %%% No need to run it. frequency distribution of genotypes in a complete landscape 
genotype = cellstr(num2str(mutant_index));
genotype = regexprep(genotype, '\W', ''); % all 1024 genotypes present

N = length(mutant_index);
ancestor_uint8 = repmat(uint8('0000000000'),N,1);
binary_genotypeID_uint8_temp = uint8(mutant_index(1,:));

for i=2:length(mutant_index) 
    binary_genotypeID_uint8_temp = [binary_genotypeID_uint8_temp; uint8(mutant_index(i,:))];
end

binary_genotypeID_uint8 = bitxor(binary_genotypeID_uint8_temp,ancestor_uint8);


% To plot the distribution of number of n-hamming distance away neighbours as a function of genotypes present

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

% Plot code

figure()

ax = gca;
initial_counts_mean = mean(mutant_counts_wFilter_wNeutral(:,17:18),2);
plot(1:N,log(initial_counts_mean));
ymax = max(initial_counts_mean);
%ylim([0 ymax+0.25*ymax]);
xlim([-25 1050]);
xticks([1 1024]);

for i=1:L
    vline(boundary(i,1)-0.5,'r:',num2str(i));
end
%vline(boundary(10,1)-0.5,'r:','10');
sub_title = sprintf('%s','Distribution of Log of Initial counts in the library');
title(sub_title);

% Plot 
figure()
bar(0:L,genotype_count_per_mutation)
ymax = max(genotype_count_per_mutation);
ylim([0  ymax+0.1*ymax]);
text(0:length(genotype_count_per_mutation)-1,genotype_count_per_mutation,num2str(genotype_count_per_mutation),'vert','bottom','horiz','center'); 
box off
title('Frequency of n-mutation genotypes')
xlabel('Number of mutations')
ylabel('Number of genotypes')

%%  To get a binary matrix with all possible combinations of placing n-1's in m total sites

 n = 2;
m = 4;
ii = nchoosek(1:m,n);
k = size(ii,1);
out = zeros(k,m);
out(sub2ind([k,m],(1:k)'*ones(1,n),ii)) = 1;


