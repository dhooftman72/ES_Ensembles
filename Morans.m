function [Morans_I,PValue,ZValue,wij] = Morans(X,Y,MaxD,VarIn,bts)
% Morans I and interated significance base on a z distribution (normal)
% Inputs are a list of (x,1) longitude values with matching (y,1)
% lattitutude values and matching variable of interest, VarIn, of length (x,1) 
% Weigthted calculation can be altered but is as default a logaritmic linear decay
% function up to maximum distance MaxD. Significance interations include
% bts runs of data permutations

%make distance matrix: X and Y need to be in the same units as MaxD, meters or kilometers
if length(X) ~= length(Y)
    display('Fatal Error, unequal grid size')
    return
end
size = length(X);
for i = 1:size
    for j = 1:size
            Dis(i,j) = sqrt (  (((X(i)-X(j))^2) + (Y(i)-Y(j))^2));  %#ok<*AGROW>
    end
end
% make weighted matrix
wij = weightedFunc(Dis,MaxD,size);
% calculate Moran's I of variable of interest
Morans_I = moransCalc(wij,VarIn,size);
% Calculate permutation distribution 
[Varsigma,Varmu] = zshuffle(wij,VarIn,size,bts);
% calculate significance
[~,PValue,~,ZValue] = ztest(Morans_I,Varmu,Varsigma);
save('all')
end

function wij = weightedFunc(Dis,MaxD,size)
% make weighted matrix, logaritmic linear function with maximum distancesz
Dist = log10(Dis+1);
wij = (log10(MaxD)-Dist)./log10(MaxD);
wij(wij<0) = 0;
for i = 1:size
    for j = 1:size
        if i == j
            wij(i,j) = 0;  %#ok<*AGROW>
        end
    end
end
end

function [Varsigma,Varmu] = zshuffle(wij,VarIn,size,bts)
if bts >= 100
    if matlabpool('size') ~= 0
        matlabpool close
    end
    matlabpool open Full
end
parfor x = 1:bts
    VarInB = VarIn((randperm(size)));        %#ok<PFBNS>
    Varbts(x) = moransCalc(wij,VarInB,size);
end
Varsigma = nanstd(Varbts);
Varmu = nanmean(Varbts);
if matlabpool('size') ~= 0
    matlabpool close
end
end

function VarOut = moransCalc(wij,VarIn,size)
MuVarIn = nanmean(VarIn);
count = 1;
for i = 1:size
    for j = 1:size
        nomi(count) = wij(i,j).*(VarIn(i)-MuVarIn).*(VarIn(j)-MuVarIn);
        count = count + 1;
    end
    denomi(i) = (VarIn(i)-MuVarIn)^2;
end
VarOut = (size./sum(sum(wij))).*(sum(nomi)./sum(denomi));
end
