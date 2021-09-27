%% Plot bar graphs of simulation results using the color palette used in Fig 3. for different hosts

% Import Color Palette and convert hex code to RGB values
 % ----- NOT USED ------

% filename = 'colorPalette_simulationResult.txt';
% color_hex = textread(filename, '%s', 'delimiter', '\n');
% total_color = length(color_hex);
% color_rgb = zeros(total_color,3);
% for i=1:total_color
%     str = color_hex{i};
%     color_rgb(i,:) = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
% end

%% 
figure3 = figure()
ax3=gca;
hold on

% Setting axes properties
ax3.YLim=[0 1];
ax3.YTick = [0:0.25:1];
ylabel('Frequency of OmpF evolution')

xlabel({'Fraction of total time \lambda population evolves on ancestral landscape' ; 'before switching to {\it malT^{-}} landscape'})
ax3.FontSize=12;
ax3.FontName='Arial';
ax3.Box='On';
linewidth_bars = 0.8;

ax3.XLim=[0 19.5];
xTicks_switch = [6:1.5:18];
xTick_malT = [3];
xTick_606 = [1]; 

hold on

% Plot host-switch simulations
Label_switch = {'0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9'};
freq_OmpF_switch = [0.333333333; 0.286666667; 0.376666667; 0.31; 0.353333333; 0.423333333; 0.406666667; 0.373333333; 0.383333333];
bar_switch = bar(xTicks_switch',freq_OmpF_switch,0.25,'BaseValue',0);
bar_switch.FaceColor = 'flat' ;%color(1);% 
% bar_switch.CData(:,:) = color_rgb(2:10,:);
bar_switch.FaceColor = [0.8 0.8 0.8];
bar_switch.LineWidth=linewidth_bars;
bar_switch.BarWidth = 0.6;

set(gca, 'XTick', xTicks_switch, 'XTickLabel', Label_switch, 'FontSize',12, 'FontName','Arial')

confidence_interval = [0.075405239; 0.057536827; 0.064622089; 0.07113716; 0.047954675; 0.081872482; 0.07771578; 0.050018177; 0.055354271];
marginOfError = confidence_interval;
h = errorbar(xTicks_switch',freq_OmpF_switch,marginOfError.','k','LineWidth',0.65);
set(h,'linestyle','none')
%%
hold on  
% plotting malT-
Label_malT = {'\it malT^{-}'};
freq_OmpF_malT = [0.186666667];
bar_malT = bar(xTick_malT,freq_OmpF_malT,0.25);
bar_malT.FaceColor=[0.8 0.8 0.8];
bar_malT.LineWidth=linewidth_bars;
bar_malT.BarWidth = 0.8;
set(gca, 'XTick', xTick_malT, 'XTickLabel', Label_malT, 'FontSize',12, 'FontName','Arial')
confidence_interval = [0.032188607];
marginOfError = confidence_interval;
h = errorbar(xTick_malT,freq_OmpF_malT,marginOfError.','k','LineWidth',0.65);
set(h,'linestyle','none')


division = xline(xTick_malT + 1.5, 'k--');
division.LineWidth = linewidth_bars + 0.75;
%%
% Plot 606 host
Label_606 = {'ancestor'};
freq_OmpF_606 = [0];
bar_606 = bar(xTick_606,freq_OmpF_606,0.25);
bar_606.FaceColor=[0.8 0.8 0.8];
bar_606.LineWidth=linewidth_bars+1.5;
bar_606.BarWidth = 0.8;
set(gca, 'XTick', xTick_606, 'XTickLabel', Label_606, 'FontSize',12, 'FontName','Arial')
confidence_interval = [0];
marginOfError = confidence_interval;
h = errorbar(xTick_606,freq_OmpF_606,marginOfError.','k','LineWidth',0.65);
set(h,'linestyle','none')

%%
% Setting global XTicks and XTick Labels and vertical partition
xTicks_all = [xTick_606 xTick_malT xTicks_switch];
malT_label = '\it malT^{-}';
Labels_all = [Label_606 Label_malT Label_switch];
set(gca, 'XTick', xTicks_all, 'XTickLabel', Labels_all, 'FontSize',7, 'FontName','Arial', 'FontWeight', 'bold')
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
print('Bar graphs for simulation result ','-dpdf','-r300')