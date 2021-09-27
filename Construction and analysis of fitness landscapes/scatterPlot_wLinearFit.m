function [p rsq_adj] = scatterPlot_wLinearFit(x, y)
%scatter(x, y)
[p, S] = polyfit(x,y,1);
%Call polyval to use the coefficients in p to predict y, naming the result yfit:
yfit = polyval(p,x);
%Compute the residual values as a vector of signed numbers:
yresid = y - yfit;
%Square the residuals and total them to obtain the residual sum of squares:
SSresid = sum(yresid.^2);
%Compute the total sum of squares of y by multiplying the variance of y by the number of observations minus 1:
SStotal = (length(y)-1) * var(y);
%Compute simple R2 for the cubic fit using the formula given in the introduction of this topic:
rsq = 1 - SSresid/SStotal;
%Finally, compute adjusted R2 to account for degrees of freedom:
rsq_adj = 1 - SSresid/SStotal * (length(y)-1)/(length(y)-length(p))

dim = [0.2 .5 .3 .3];
str = {sprintf('%s', 'Linear fit - '),sprintf('%s %.4f', 'Slope -', p(1)), sprintf('%s %.4f', 'R^{2} -', rsq_adj)};
annotation('textbox',dim,'String',str,'FitBoxToText','on');

% to plot the line fitted
hold on;
plot(x,yfit,'r-.');

end