% Plot bar graphs for difference betweennn same host switching results and two host switching results

% DATA from 'Final simulation analysis_fitness landscape.xlsx'

%%  For figure ecc4 to ecc4 vs 606 to ecc4 switching. 
% Sheet 'Ecc4 to Ecc4, N divided by 10' and 'N divided by 10, T'

fig_1 = figure()
ax1=gca;
hold on
xTicks = [1;2];
linewidth_bars = 0.8;

% Plot host-switch simulations
Label_xaxis = {'\it malT^{-} to \it malT^{-}', 'ancestor to \it malT^{-}'};
freq_OmpF = [0.192883895; 0.305925926];
bar_ecc4 = bar(xTicks,freq_OmpF,0.25);%,'BaseValue',0);
bar_ecc4.LineWidth=linewidth_bars;
% bar_1v0.BarWidth = 0.6;

stdError=[0.052767985; 0.066641922];
marginOfError = stdError;
h = errorbar(xTicks,freq_OmpF,marginOfError.','k','LineWidth',0.65);
set(h,'linestyle','none')

% Setting axes properties
ax1.YLim=[0 0.75];
ax1.YTick = [0:0.25:0.75];
ax1.XTick = xTicks;
ax1.XTickLabel = Label_xaxis;
ylabel('Frequency of OmpF evolution')
ax1.FontSize=12;
ax1.FontName='Helvetica';

%%  For figure 606 to 606 vs 606 to ecc4 switching.
% Sheet '606 to 606, common, N by 10' and 'Common genotypes, N over 10, T'
fig_2 = figure()
ax1=gca;
hold on
xTicks = [1;2];
linewidth_bars = 0.8;

% Plot host-switch simulations
Label_xaxis = {'\it ancestor to ancestor', 'ancestor to \it malT^{-}'};
freq_OmpF = [0.212962963; 0.337777778];
bar_606 = bar(xTicks,freq_OmpF,0.25);%,'BaseValue',0);
bar_606.LineWidth=linewidth_bars;
% bar_1v0.BarWidth = 0.6;

stdError=[0.057935793; 0.080502828];
marginOfError = stdError;
h = errorbar(xTicks,freq_OmpF,marginOfError.','k','LineWidth',0.65);
set(h,'linestyle','none')

% Setting axes properties
ax1.YLim=[0 0.75];
ax1.YTick = [0:0.25:0.75];
ax1.XTick = xTicks;
ax1.XTickLabel = Label_xaxis;
ylabel('Frequency of OmpF evolution')
ax1.FontSize=12;
ax1.FontName='Helvetica';

