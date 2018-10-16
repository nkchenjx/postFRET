function [binData, DataBinLength] =binning(trajData, Datalength, binsize)

numMol = size(Datalength,2);

DataRange = Datalength;
DataBinLength = zeros(1, numMol);
binData = [];
DataBinLength = [];


if numMol == 1
    DataBinLength = floor(size(trajData,2)/binsize);
    DataMolChop = trajData(1:DataBinLength*binsize);
    bintemp = DataMolChop; 
    bintemp = reshape(bintemp, binsize, DataBinLength);
    binData = sum(bintemp,1);  

else
    for i = 2:numMol
      DataRange(i) = sum(Datalength(1:i));
    end

    for i = 1:numMol % bin each molecule
        if i == 1
            DataMol = trajData(1:Datalength(i));
        else
            DataMol = trajData(DataRange(i-1)+1:DataRange(i));
        end
        DataBinLength(i) = floor(size(DataMol,2)/binsize); % new length of the data
        DataMolChop = DataMol(:, 1:DataBinLength(i)*binsize); % choped into multi of binsize

        bintemp = DataMolChop; 
        bintemp = reshape(bintemp, binsize, DataBinLength(i));
        DataBinMol = sum(bintemp,1);  
        binData = [binData,DataBinMol];
    end
end