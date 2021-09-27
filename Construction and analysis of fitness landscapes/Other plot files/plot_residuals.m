%% Plot std dev for 606
figure()
subplot(2,1,1)
selectionRate_wGenotypeID_606 = selectionRate_wGenotypeID_WT;
error_606 = nanstd(selectionRate_wGenotypeID_606(:,2:5),0,2);
mean_606 = nanmean(selectionRate_wGenotypeID_606(:,2:5),2);
scatter(mean_606,error_606, 10, 'filled', 'k')
non_zero_indices_606 = find(error_606);
ax = gca;
hold on
scatterPlot_wLinearFit(mean_606(non_zero_indices_606),error_606(non_zero_indices_606))
xlabel('Average selection rate')
ylabel('Standard deviation')
title('ancestor')
ax.FontSize=12;
ax.FontName='Arial';

% Plot std dev for malT
subplot(2,1,2)
error_malT = nanstd(selectionRate_wGenotypeID_malT(:,2:5),0,2);
mean_malT = nanmean(selectionRate_wGenotypeID_malT(:,2:5),2);
scatter(mean_malT,error_malT, 10, 'filled', 'k')
non_zero_indices_malT = find(error_malT);
ax = gca;
hold on
scatterPlot_wLinearFit(mean_malT(non_zero_indices_malT),error_malT(non_zero_indices_malT))
xlabel('Average selection rate')
ylabel('Standard deviation')
title('\it malT^{-}')
ax.FontSize=12;
ax.FontName='Arial';
ax.XLim = [-1 11];
ax.XTick = [-1:1:11]

%% Plot with colorbar according to the 
figure()
x = phageInitial(selectionRate_wGenotypeID_WT(:,1)+1);
y = x./max(x);
y = log(y);
pointsize = 20;
%scatter(commonID_WT_fitness_bc1, commonID_WT_fitness_bc2, pointsize,y,'filled')%'MarkerEdgeColor', 'black');
c = colorbar('Direction','normal');
scatter(mean_606,error_606, pointsize,y, 'filled')
%%
z = x(non_zero_indices_606);
mean_606_non_zero = mean_606(non_zero_indices_606);
error_606_non_zero = error_606(non_zero_indices_606);
[z_sort,sortIdx] = sort(z,'descend');
% sort B using the sorting index
mean_606_sort = mean_606_non_zero(sortIdx);
error_606_sort = error_606_non_zero(sortIdx);

%% Plot residuals

residual_606 = selectionRate_wGenotypeID_606(:,2:5)-mean_606;
residual_mean_606 = nanmean(residual_606,2);
figure()
scatter(mean_606,residual_606(:,1))
hold on
scatter(mean_606,residual_606(:,2))
scatter(mean_606,residual_606(:,3))
scatter(mean_606,residual_606(:,4))
xlabel('Mean selection rate across 4 replicates')
title('606')
ylabel('Residuals- r_{i} = x_{i} - x_{average}')