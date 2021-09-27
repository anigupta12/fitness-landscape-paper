%% Last edited- 30Jan'17

function [peakStruct] = localPeakAnalysis(n, genotypeID, genotype, selectionRates, localPeak)

% Returns a structure, named peakStruct, which includes-
% a. genotype (cell) -  binary representation of genotypes which are peaks detected in simple nearest neighbour interaction;
% b. selectionCoeff (2xno.ofPeaks cell)- contains genotypes in 1st row, and all the selection coefficient of it (1st element) and its neighbours in the second row as  column vector 
% c. totalNeighbours (scalar): total no. of neighbours above detection threshold for the peaks measured

% Input- 
% a. n (scalar total no. of point mutation sites)
% b. genotype (totalNo.Genotypesx1 cell)- binary representation of genotypes for which data is present, i.e. above the detection threshold
% c. genotypeID (double totalGenotypesx1 matrix) - decimal representation of above genotypes
% d. selectionRates (double totalGenotypesx4 matrix) - selection Rates of the above genotypes
% e. localPeaks (binary totalGenotypesx1)- genotypes which are local peaks. peaks have value of 1 and others 0

peakStruct = struct; % creating a structure
a= find(localPeak==1); % 'a' is a nx1 matrix where n is the the number of peaks, and the entries of 'a' are the index no. of the local peaks in the imported data (also b=a as local peaks is a column vector)
for u=1:length(a)
    tempGenotype = genotype{a(u),1};
    peakStruct(u).genotype = genotype(a(u),1); % this creates a cell
    peakStruct(u).selectionCoeff{1,1} = genotype{a(u),1}; 
    %genotype_decimal = bin2dec(tempGenotype);
    peakStruct(u).selectionCoeff{2,1} = transp(selectionRates(a(u),:)); % 1st row of the selection coeff is the genotype peak's selection rate
    %peakStruct(u).selectionCoeff{2,1} = transp(selectionRates(genotype_decimal+1,:)); 
    
    neighboursPresent = 0;
    for j = 1:n
        neighbour = tempGenotype;
        newBit = mod(bin2dec(neighbour(j))+1,2);
        neighbour(j) = num2str(newBit);
        %disp(neighbour)
        neighbourID = bin2dec(neighbour);
        neighbourID_index = find(genotypeID==neighbourID);
        %disp(neighbourID_index)
        if isempty(neighbourID_index)==false
            neighboursPresent = neighboursPresent +1;
            peakStruct(u).selectionCoeff{1,1+neighboursPresent} = neighbour; % adding binary representation of the neighbouring genotype
            peakStruct(u).selectionCoeff{2,1+neighboursPresent} = transp(selectionRates(neighbourID_index,:)); % adding selection rate of the neighbouring genotype
        %peakStruct(u).selectionCoeff{2,1+neighboursPresent} = transp(selectionRates(neighbourID+1,:)); 
        end
        %adder = dec2bin(2^j,n);
        %neighbour = dec2bin(mod(bin2dec(genotype{i,1}) + bin2dec(adder),1024),10);  %to add numbers in binary
    end
    peakStruct(u).totalNeighbours = neighboursPresent;
end
end
