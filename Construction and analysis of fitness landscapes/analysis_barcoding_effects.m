%% To test barcoding effects 

% Using mutant_counts_wNeutral_wFilter
data = mutant_counts_wNeutral_wFilter;
NaNindices=data(:)==0.1;  % Removing all the entries which are 0.1 and represent 0 actually
data_NaN = data;
data_NaN(NaNindices)=nan;   % Replacing with NaN
%% TO SEE BARCODING EFFECTS IN RAW COUNTS

start=17; % 1-8 for 606, 9-16 for ECC4, and 17-18 for Initial
last=18;
p = zeros(ceil((last-start)/2),1);
rsq_adj = zeros(ceil((last-start)/2),1);
count=0;
for i = start:2:last
    count=count+1;
    genotypeID = [0:1023]';
    
    bc1 = data_NaN(:,i);    
    bc2 = data_NaN(:,i+1);
    fitness_both_bc = [bc1 bc2];
    fitness_sum = sum(fitness_both_bc,2);
    [a, b] = find(isnan(fitness_sum));
    bc1(a) = [];
    bc2(a) = [];
    genotypeID_filtered=genotypeID;
    genotypeID_filtered(a) = [];
    
    
    bc1_filtered = [genotypeID_filtered ,bc1];
    bc2_filtered = [genotypeID_filtered ,bc2];
    
    l = length(bc1_filtered);
    bc1_filtered(l-1:l,:)=[];
    bc2_filtered(l-1:l,:)=[];
    
    [commonID, commonID_fitness_bc1, commonID_fitness_bc2] = compare_fitness_common_mutants(bc1_filtered,bc2_filtered);
    [p_temp, rsq_adj_temp] = scatterPlot_wLinearFit(commonID_fitness_bc1, commonID_fitness_bc2);
    p(count,1) = p_temp(1);
    rsq_adj(count,1) = rsq_adj_temp;
end

title('Initial')
xlabel('Barcode 1')
ylabel('Barcode 2') 

figure()
x = ([1])';
scatter( x, p, 'filled');
hold on
scatter( x, rsq_adj, 'filled');
legend('Slopes','R^{2}')
title('Slopes for Initial')
xlabel('Replicates')
ylabel('Slope of fitted line')

% dim = [0.2 .5 .3 .3];
% str = {sprintf('%s', 'Linear fit - '),sprintf('%s %.4f', 'Slope -', p_WT(1)), sprintf('%s %.4f', 'R^{2} -', rsq_adj_WT)};
% annotation('textbox',dim,'String',str,'FitBoxToText','on');

%% TO SEE BARCODE EFFECTS IN SELECTION RATE

dataZero=data;
dataZero(NaNindices)=0;

phageInitial_bc1=dataZero(:,17);
phageInitial_bc2=dataZero(:,18);
%phageInitial_NaN=mean(dataNaN(:,17:18),2);    % mean of 1 and NaN will be NaN

phageWTa_bc1=dataZero(:,1);
phageWTb_bc1=dataZero(:,3);
phageWTc_bc1=dataZero(:,5);
phageWTd_bc1=dataZero(:,7);
phageECCa_bc1=dataZero(:,9);
phageECCb_bc1=dataZero(:,11);
phageECCc_bc1=dataZero(:,13);
phageECCd_bc1=dataZero(:,15);


phageWTa_bc2=dataZero(:,2);
phageWTb_bc2=dataZero(:,4);
phageWTc_bc2=dataZero(:,6);
phageWTd_bc2=dataZero(:,8);
phageECCa_bc2=dataZero(:,10);
phageECCb_bc2=dataZero(:,12);
phageECCc_bc2=dataZero(:,14);
phageECCd_bc2=dataZero(:,16);

phageAll_bc1=[phageWTa_bc1 phageWTb_bc1 phageWTc_bc1 phageWTd_bc1 phageECCa_bc1 phageECCb_bc1 phageECCc_bc1 phageECCd_bc1];
phageAll_bc2=[phageWTa_bc2 phageWTb_bc2 phageWTc_bc2 phageWTd_bc2 phageECCa_bc2 phageECCb_bc2 phageECCc_bc2 phageECCd_bc2];

[ ~, ~, selectionRate_wGenotypeID_WT_bc1, selectionRate_mean_wGenotypeID_WT_bc1] = selectionRate(phageAll_bc1, phageInitial_bc1, 1, 4, 1);
[ ~, ~, selectionRate_wGenotypeID_malT_bc1, selectionRate_mean_wGenotypeID_malT_bc1] = selectionRate(phageAll_bc1, phageInitial_bc1, 5, 8, 1);

[ ~, ~, selectionRate_wGenotypeID_WT_bc2, selectionRate_mean_wGenotypeID_WT_bc2] = selectionRate(phageAll_bc2, phageInitial_bc1, 1, 4, 1);
[ ~, ~, selectionRate_wGenotypeID_malT_bc2, selectionRate_mean_wGenotypeID_malT_bc2] = selectionRate(phageAll_bc2, phageInitial_bc1, 5, 8, 1);


%%%%%%%%%--------%%%%%%%%%%%%%%%
%% To get the number of filtered mutants and their identity

[count_total_mutants_WT, count_filtered_mutants_WT, filtered_mutants_ID_WT, filtered_mutant_cell_WT] = compare_unique_mutants(selectionRate_mean_wGenotypeID_WT_bc1(:,1),selectionRate_mean_wGenotypeID_WT_bc2(:,1));
[count_total_mutants_malT, count_filtered_mutants_malT, filtered_mutants_ID_malT, filtered_mutant_cell_malT] = compare_unique_mutants(selectionRate_mean_wGenotypeID_malT_bc1(:,1),selectionRate_mean_wGenotypeID_malT_bc2(:,1));

%% To plot the above data
combined_filtered_mutant_cell = { filtered_mutant_cell_WT{:,:}; filtered_mutant_cell_malT{:,:}};
combined_filtered_mutant_cell(:,4) = {'606', 'malT'}; 
plot_filtered_mutants_table(combined_filtered_mutant_cell)


% Results- Neutral filtering is stronger than just filtering. So we are
% good

%% To plot difference in fitness due to filtering

[commonID_WT, commonID_WT_fitness_bc1, commonID_WT_fitness_bc2] = compare_fitness_common_mutants(selectionRate_mean_wGenotypeID_WT_bc1, selectionRate_mean_wGenotypeID_WT_bc2);
[commonID_malT, commonID_malT_fitness_bc1, commonID_malT_fitness_bc2] = compare_fitness_common_mutants(selectionRate_mean_wGenotypeID_malT_bc1, selectionRate_mean_wGenotypeID_malT_bc2);
%% Plot for 606
figure()
%scatter(commonID_WT_fitness_wFilter, commonID_WT_fitness_woFilter)

[p_WT, rsq_adj_WT] = scatterPlot_wLinearFit(commonID_WT_fitness_bc1, commonID_WT_fitness_bc2);
% Turn off the scatter plot in above scatterPlot_wLinearFit so that you can
% replot it with different settings as below

hold on
ax606 = gca;
% To log normalize for colorbar according to initial counts 
x = phageInitial(commonID_WT+1);
y = x./max(x);
y = log(y);
pointsize = 20;
scatter(commonID_WT_fitness_bc1, commonID_WT_fitness_bc2, pointsize,y,'filled')%'MarkerEdgeColor', 'black');
c = colorbar('Direction','normal');
title('REL606')
ax606.XGrid = 'on';
ax606.YGrid = 'on';
box(ax606,'on')
xlabel('Selection rate using ''Barcode 1''')
ylabel('Selection rate using ''Barcode 2''')
title('REL606')
%dim = [0.2 .5 .3 .3];
% After this, all fonts bold, size-14, and style- Helvetica

%print('5by3DimensionsFigure','-dtiff','-r300')
%% Plot for malT
figure()
%scatter(commonID_malT_fitness_wFilter, commonID_malT_fitness_woFilter)

[p_malT, rsq_adj_malT] = scatterPlot_wLinearFit(commonID_malT_fitness_bc1, commonID_malT_fitness_bc2);
% Turn off the scatter plot in above scatterPlot_wLinearFit so that you can
% replot it with different settings as below

hold on
axmalT = gca;
% To log normalize for colorbar according to initial counts 
x = phageInitial(commonID_malT+1);
y = x./max(x);
y = log(y);
pointsize = 20;
scatter(commonID_malT_fitness_bc1, commonID_malT_fitness_bc2, pointsize,y,'filled')%'MarkerEdgeColor', 'black');
c = colorbar('Direction','normal');
title('RELmalT')
axmalT.XGrid = 'on';
axmalT.YGrid = 'on';
box(axmalT,'on')
xlabel('Selection rate using ''Barcode 1''')
ylabel('Selection rate using ''Barcode 2''')
title('EcC4')
%dim = [0.2 .5 .3 .3];
% After this, all fonts bold, size-14, and style- Helvetica

%print('5by3DimensionsFigure','-dtiff','-r300')
%%
figure()
[p_malT, rsq_adj_malT] = scatterPlot_wLinearFit(commonID_malT_fitness_bc1, commonID_malT_fitness_bc2);
xlabel('Selection rate using ''Barcode 1''')
ylabel('Selection rate using ''Barcode 2''')
title('EcC4')
dim = [0.2 .5 .3 .3];
str = {sprintf('%s', 'Linear fit - '),sprintf('%s %.4f', 'Slope -', p_malT(1)), sprintf('%s %.4f', 'R^{2} -', rsq_adj_malT)};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
%% To plot 1:1 line on the plots

hold on
xl=xlim;
yl=xl;
plot(xl,yl,'-k','LineWidth',0.5)
%% Rough space to figure out final plot settings

figure()

hold on
%xl = xlim;
%yl = xl;
%plot(xl,yl)
ax606 = gca;
% set(gca,'FontSize',14)
% set(gca,'FontName','Helvetica')
% set(gca,'FontWeight','bold')
x = phageInitial(commonID_WT+1);
%y = max(y) - y;
y = x./max(x);
y = log(y);
pointsize = 20;
scatter(commonID_WT_fitness_bc1, commonID_WT_fitness_bc2, pointsize,y,'filled')%'MarkerEdgeColor', 'black');
c = colorbar('Direction','normal');
hold on

xlabel('Selection rate using ''Barcode 1''')
ylabel('Selection rate using ''Barcode 2''')
title('REL606')
ax606.XGrid = 'on';
ax606.YGrid = 'on';
%xticks(ax606,(-3:1:7))
%xlim(ax606, [-3.25 7.25])
%ylim(ax606, [-3.25 7.25])
box(ax606,'on')
%dim = [0.2 .5 .3 .3];
