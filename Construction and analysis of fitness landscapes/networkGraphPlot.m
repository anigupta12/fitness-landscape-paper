%% Last edit- 2/17/17

% Output- Plots two network graphs- full, and one with only peaks and paths to them
% Input-
% a. graphData (structure) - output of networkGraph.m file. Check that for
% more info
% b. l (scalar)- total no. of genotypes present in the data
% c. treatment (string)- denoting which treatment it is for title of the
% graph
function [u] = networkGraphPlot(graphData, l, treatment)
% For full network graph
figure('Name',treatment)
box on
for u = 1:l
    hold on
    plot(graphData(u).xdata, graphData(u).ydata,'-o','LineWidth', 0.2,'Color',[0 0 0.9],'MarkerSize', 3, 'MarkerFaceColor',[0 0 0],'MarkerEdgeColor' ,[0 0 0.5])

    %To plot ANOVA local peaks, uncomment below if, else loops. NOTE these
    %are not bonferroni passed peaks

%     if graphData(u).localPeak==1
%         plot(graphData(u).xdata, graphData(u).ydata,'r-','LineWidth', 0.75)
%         plot(graphData(u).xdata(1,1), graphData(u).ydata(1,1),'ro', 'MarkerSize', 8,'MarkerFaceColor','r')
%     else
%     end

end
hold on
ylabel('Selection Rate');
xlabel('Number of Mutations');
title(treatment,'FontSize',14, 'FontName', 'Helvetica');
hold on
%ylim([-3 11]);
plot(xlim,[0 0],'--k','LineWidth',0.75);

%% TO PLOT Network Graph with only Local Peaks
% figure('Name',strcat(treatment,' local peaks'))
% box on
% for u=1:l
%     hold on
%     if graphData(u).localPeak==1
%         plot(graphData(u).xdata, graphData(u).ydata,'r-','LineWidth', 0.75)
%         plot(graphData(u).xdata(1,1), graphData(u).ydata(1,1),'ro', 'MarkerSize', 8,'MarkerFaceColor','r')
%     end
%     
% end
% hold on
% ax=gca;
% ax.XLim=[0 10];
% ylabel('Fitness');
% xlabel('# of Mutations');
% title(strcat(treatment,' local peaks'),'FontSize',14);
% plot(xlim,[0 0],'--k','LineWidth',0.35);
% end
