%% To edit properties of already drawn MATLAB figures

fig1 = get(gcf);
%%
set(findall(gcf,'-property','FontName'),'FontName','Helvetica')
%%
set(findall(gcf,'-property','FontWeight'),'FontWeight','normal')
%fig1.Children(1).YTick = [0:0.2:2]