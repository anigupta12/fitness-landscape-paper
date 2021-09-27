%%  Last edited- 30Jan'17
function [anovaOutput,  validAnovaPeaks, multiCompareOutput, validPeaks, finalpassedGenotype] = anova(peakStruct, totalLocalPeaks, test, alpha )

%Output- 
% a. anovaOutput (1xtotalLocalPeaks cell) - output of anova 
% b. validAnovaPeaks (scalar)- No. of peaks which passed anova
% c. multiCompare (validAnovaPeaksx1 cell) - output of tukey's or bonferrni test
% d. validPeaks (scalar)- No. of peaks which passed above multicompare test
% e. finalPassedGenotype (validPeaksx2 cell)- 1st column are genotypes which passed multicompare and 2nd column are the corresponding index no. in peakStruct

%Input- 
% a. peakStruct (structutre)- output of localPeakAnalysis.m, input for anova function
% b. totalLocalPeaks- total no. of local peaks
% c. test (string)- can be tukey or bonferroni or some other valid input of anova function
% d. alpha(scalar)- p-value or alpha value for the test

validAnovaPeaks = 0;
validPeaks = 0;
t=1;    % just a counter
anovaOutput = cell(totalLocalPeaks,1);
multiCompareOutput = cell(totalLocalPeaks,1);
finalpassedGenotype = cell(totalLocalPeaks,2);

for u=1:totalLocalPeaks
    [p, table, stats] = anova1(cell2mat(peakStruct(u).selectionCoeff(2,:)),[],'off');
    anovaOutput{u} = {p, table, stats};
    if anovaOutput{u}{1} <= 0.05  % checking if the p-value after ANOVA is less than 0.05
        
        if anovaOutput{u}{3}.df ~= 0
            validAnovaPeaks = validAnovaPeaks + 1;
            peakStruct(u).validAnovaPeak = 1;
            multiCompareOutput{u} = multcompare(anovaOutput{u}{3}, 'Alpha',alpha,'CType',test);   % default is tukey's test        
            disp(u)
        % to check if the local peak is significantly different (greater) than others
            significantPValues = multiCompareOutput{u}(1:peakStruct(u).totalNeighbours,6) > alpha;
            if sum(significantPValues) == 0
                validPeaks = validPeaks + 1;
                peakStruct(u).validPeak = 1;
                disp(peakStruct(u).genotype)
                finalpassedGenotype{t,1} = [peakStruct(u).genotype];
                finalpassedGenotype{t,2} = u; 
                t=t+1;
            end      
        end
           
    end
    % multicompare{i} = {multcompare(anova{i}{3}, 'Alpha',0.05)};   % default is tukey's test
end
multiCompareOutput(~cellfun('isempty', multiCompareOutput));    % concatenating empty cells if any, since it was initialized earlier with maximum no. of possible statistically valid peaks
finalpassedGenotype(~cellfun('isempty', finalpassedGenotype));  
end