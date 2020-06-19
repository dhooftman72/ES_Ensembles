function [Points,Regressions,AutoCorrelation] = PointRegression(Deviations,Setnames)
clc
warning off
close all
%load('Ensembles.mat')
btsfact = 5;
Service_nots = unique(Deviations.Service,'stable');
whichones = (1:length(Deviations.Service));%
Points = dataset(Deviations.Service(whichones),'Varnames','Services');
Points.Dataset = Deviations.SetID;
Points.SEM = Deviations.MeanStd (whichones);
CorFac = Deviations.N(whichones);
Points.SEM = Points.SEM./sqrt(CorFac);
Points.Means = Deviations.Mean(whichones);
Points.SEMLogit = log10((Points.SEM./(1-Points.SEM))+1);
Points.SEMLogit(Points.SEMLogit>0.99) = 0.99;
Size = length(Points.SEM);
AutoCorrelation = dataset({'OverallAccuracy'},'Varnames', {'Set'});

str = sprintf('      MoransI for set = %s ',char('OverallAccuracy'));
disp(str)
[Morans_I,PValue,ZValue,wij] = Morans(Deviations.Long,Deviations.Lat,1000000,Points.Means,(Size.*btsfact)); 
AutoCorrelation = Autoset(AutoCorrelation,Size,{'OverallAccuracy'},Morans_I,PValue,ZValue,1); 

str = sprintf('      MoransI for set = %s ',char('SEM'));
 disp(str)
[Morans_I,PValue,ZValue,~] = Morans(Deviations.Long,Deviations.Lat,1000000,Points.SEM,(Size.*btsfact)); 
AutoCorrelation = Autoset(AutoCorrelation,Size,{'SEMUncertainty'},Morans_I,PValue,ZValue,2); 

str = sprintf('      MoransI for set = %s ',char('OverallUncertaintyLogit'));
disp(str)
[Morans_I,PValue,ZValue,~] = Morans(Deviations.Long,Deviations.Lat,1000000,Points.SEMLogit,(Size.*btsfact)); 
AutoCorrelation = Autoset(AutoCorrelation,Size,{'OverallUncertaintyLogit'},Morans_I,PValue,ZValue,3); 

for i = 1:1:16
    display(' ')
    display(' ')
    str = sprintf('      MoransI for set = %s ',char(Setnames(i)));
    disp(str)
    test = find(Points.Dataset == i);
    [Morans_I,PValue,ZValue,~] = Morans(Deviations.Long(test),Deviations.Lat(test),1000000,Points.Means(test),(length(Points.Means(test)).*btsfact));
    AutoCorrelation = Autoset(AutoCorrelation,length(Points.Means(test)),Setnames(i),Morans_I,PValue,ZValue,(i+3));
end
if exist('wij','var') ~= 1
    Auto = Deviations.Autowij;
else
    Y = Points.Means;
    for i = 1:Size
        Autot = 0;
        for j = 1:Size
            if i ~= j
                Autot = Autot + (wij(i,j).*Y(j));
            end
        end
        Auto(i,1) = Autot./sum(wij(i,:)); %#ok<*AGROW>
    end
    
end
Points.Auto(:,1) = Auto;
% Interaction model
[~,outs,stats] = anovan(Points.Means,{Auto,Points.Services,Points.SEMLogit},'sstype',1,...
     'model',[1 0 0 ; 0 1 0 ; 0 0 1  ; 0 1 1],'continuous', [1,3],'display', 'on',...
    'varnames', {'AutoCorrelation','Service','Uncertainty'});
Points.MeanCorrected =  Points.Means- stats.resid;
Regressions.Mean = outs;
Regressions.Stats = stats;

Regressions.EffectPerES = dataset(stats.coeffs((length(stats.coeffs))-(5)),'Varnames',char(Service_nots(1)));
for i = 5:-1:1
    Regressions.EffectPerES.(genvarname(char(Service_nots(7-i)))) = stats.coeffs((length(stats.coeffs))-(i-1));
end
[Morans_I,PValue,ZValue,~] = Morans(Deviations.Long,Deviations.Lat,1000000,Points.MeanCorrected,(Size.*btsfact));
AutoCorrelation = Autoset(AutoCorrelation,Size,{'MarginalValues'},Morans_I,PValue,ZValue,20); 
[Morans_I,PValue,ZValue,~] = Morans(Deviations.Long,Deviations.Lat,1000000,stats.resid,(Size.*btsfact));
AutoCorrelation = Autoset(AutoCorrelation,Size,{'Residuals'},Morans_I,PValue,ZValue,21); 

expect = Points.MeanCorrected;
meanY = mean(Points.Means);
for t = 1:1:Size
    ssres(t) = ((expect(t)-Points.Means(t)).^2);
    sstot(t) =  ((Points.Means(t)-meanY).^2);
end
Regressions.Mean(1,8) = {'R2'};
Regressions.Mean(2,8) = {1- ((sum(ssres)) /(sum(sstot)))};
end

function AutoCorrelation = Autoset(AutoCorrelation,size,txt,Morans_I,PValue,ZValue,i)
    AutoCorrelation.Set(i,1) = txt; 
    AutoCorrelation.N(i,1) = size;
AutoCorrelation.MoransI(i,1) = Morans_I;
AutoCorrelation.ZValue(i,1) = ZValue;
AutoCorrelation.PValue(i,1) = PValue;
end
