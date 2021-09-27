%

figure1 = figure()

ax3=gca;
hold on
axis off
ax3.YLim=[0 1];
ax3.YTick = [];
 
xTicks_switch = [6:1.5:18];
xTick_malT = [3];
xTick_606 = [1]; 
xTicks_all = [xTick_606 xTick_malT xTicks_switch];



ax3.XLim=[0 19.5];

baseLine_cartoon = 0.6;
ratio_time = [ 1 0; 0 1; 0.1 0.9; 0.2 0.8; 0.3 0.7; 0.4 0.6; 0.5 0.5; 0.6 0.4; 0.7 0.3; 0.8 0.2; 0.9 0.1];
y_time = ratio_time.*0.2;
y_time(:,1) = y_time(:,1) + baseLine_cartoon;
baseLine_cartoon = 0.6;
ratio_time = [ 1 0; 0 1; 0.1 0.9; 0.2 0.8; 0.3 0.7; 0.4 0.6; 0.5 0.5; 0.6 0.4; 0.7 0.3; 0.8 0.2; 0.9 0.1];
y_time = ratio_time.*0.2;
y_time(:,1) = y_time(:,1) + baseLine_cartoon;
%%
xTicks_all = xTicks_all';
bar_time = bar(xTicks_all(3:end), y_time(3:end,:),'stacked' , 'BaseValue', baseLine_cartoon, 'ShowBaseLine', 'off','EdgeColor',[1 1 1], 'BarWidth', 0.6);
bar_time(2).FaceColor = '#F9ED32';
bar_time(1).FaceColor = '#27AAE1';

%%

bar_time = bar(xTicks_all(1:2), y_time(1:2,:),'stacked' , 'BaseValue', baseLine_cartoon, 'ShowBaseLine', 'off','EdgeColor',[1 1 1], 'BarWidth', 0.6);
bar_time(2).FaceColor = '#F9ED32';
bar_time(1).FaceColor = '#27AAE1';
%%
ax3 = gca; 
ax3.XTickMode = 'manual';
ax3.YTickMode = 'manual';
ax3.ZTickMode = 'manual';
ax3.XLimMode = 'manual';
ax3.YLimMode = 'manual';
ax3.ZLimMode = 'manual';

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 3];
print('cartoon1','-dpdf','-r300')

%set(gca, 'FontSize',7, 'FontName','Arial', 'FontWeight', 'bold')
