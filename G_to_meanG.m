% Function: G_to_meanG belonging to Willcock et al.
%This is an Example for one ES for 6 models to be translated to mean and 
%SEM grid maps and perctile maps of hot and cold areas 
%%
% Mean, Median and SEM variation among maps of ES
% Load and winzorese models: do this for all
%Make sure all maps are same size and oevrlay in ArcGIS
load('Model_1_grid.asc','-ascii')
perc95 = prctile(Model_1_grid(Model_1_grid>-1),95);
Model_1_grid(Model_1_grid== -9999) = nan;
Invest = Model_1_grid./perc95;
Invest(Model_1>1) = 1;
clear Model_1_grid
.......
    
x_max = size(Model_1,1);
y_max = size(Model_1,2);
% Translate models to 5 percentiles
% Do this for all models
Percentiles.Model_1 = zeros(x_max, y_max);
Percentiles.Model_1 = Percentiles.Model_1 + 3;  
A = reshape(Model_1,1,(x_max* y_max));
percs = prctile(A(A>-1),[10,25,75,90]);
clear A
Percentiles.Model_1(Model_1<= percs(2)) = 2;
Percentiles.Model_1(Model_1<= percs(1)) = 1;
Percentiles.Model_1(Model_1>= percs(3)) = 4;
Percentiles.Model_1(Model_1>= percs(4)) = 5;
Percentiles.Model_1(isnan(Model_1)==1) = -9999;
.......
    
% Count Models per gridcell
Percentiles.zeros = zeros(x_max, y_max);
Perc_tmp1_zero = Percentiles.zeros;
Perc_tmp2_zero = Percentiles.zeros;
Perc_tmp3_zero = Percentiles.zeros;
Perc_tmp4_zero = Percentiles.zeros;
Perc_tmp5_zero = Percentiles.zeros;
Perc_tmp6_zero = Percentiles.zeros;
Perc_tmp1_zero(Percentiles.Model_1 >-1) = 1; 
Perc_tmp2_zero(Percentiles.Model_2 > -1) = 1; 
Perc_tmp3_zero(Percentiles.Model_3> -1) = 1; 
Perc_tmp4_zero(Percentiles.Model_4 > -1 ) = 1; 
Perc_tmp5_zero(Percentiles.Model_5 > -1) = 1; 
Perc_tmp6_zero(Percentiles.Model_6 > -1) = 1; 
Percentiles.Amount = bsxfun(@plus,Perc_tmp1_zero,Perc_tmp2_zero);
Percentiles.Amount = bsxfun(@plus,Percentiles.Amount,Perc_tmp3_zero);
Percentiles.Amount = bsxfun(@plus,Percentiles.Amount,Perc_tmp4_zero);
Percentiles.Amount = bsxfun(@plus,Percentiles.Amount,Perc_tmp5_zero);
Percentiles.Amount = bsxfun(@plus,Percentiles.Amount,Perc_tmp6_zero);
Percentiles.Outside = zeros(x_max, y_max);
Percentiles.Outside(Percentiles.Amount <=1) = -9999;
Percentiles = rmfield(Percentiles,'zeros');
clear Perc_tmp*

% Mena, Median and Variation per grid cell for Figure 2
display('Sem preparation')
Afull(:,:,1) = Model_1;
Afull(:,:,2) = Model_2;
Afull(:,:,3) = Model_3;
Afull(:,:,4) = Model_4;
Afull(:,:,5) = Model_5;
Afull(:,:,6) = Model_6;
Afulltmp = Afull;
Afulltmp(isnan(Afull)==1) = -9999;
Model_count = size(Afulltmp,3)-sum( Afulltmp<0,3);
clear Afulltmp
means = nanmean(Afull,3);
means(Percentiles.Outside == -9999) = -9999;
save('Mean.asc', 'means', '-ascii')
medians = nanmedian(Afull,3);
medians(Percentiles.Outside == -9999) = -9999;
save('Median.asc', 'medians', '-ascii')
stdes = nanstd(Afull,0,3);
stdes(Percentiles.Outside == -9999) = -9999;
save('Std.asc', 'stdes', '-ascii')
sem= (nanstd(Afull,0,3))./sqrt(Model_count);
sem(Percentiles.Outside == -9999) = -9999;
save('SEM_grazing.asc', 'sem', '-ascii')
clear stdes sem Afull medians
%clear means !! is used later

%% Count percentiles (SI 2 of Willcock et al)
%10%
display('combined 10%')
Perc_tmp1 = Percentiles.Outside;
Perc_tmp2 = Percentiles.Outside;
Perc_tmp3 = Percentiles.Outside;
Perc_tmp4 = Percentiles.Outside;
Perc_tmp5 = Percentiles.Outside;
Perc_tmp6 = Percentiles.Outside;
Perc_tmp1(Percentiles.Model_1 == 1) = 1; 
Perc_tmp2(Percentiles.CN == 1) = 1; 
Perc_tmp3(Percentiles.LPJ == 1) = 1; 
Perc_tmp4(Percentiles.Benefit == 1) = 1; 
Perc_tmp5(Percentiles.ScholesInt == 1) = 1; 
Perc_tmp6(Percentiles.ScholesSA == 1) = 1; 
Percentiles.Cool10 = bsxfun(@plus,Perc_tmp1,Perc_tmp2);
Percentiles.Cool10 = bsxfun(@plus,Percentiles.Cool10,Perc_tmp3);
Percentiles.Cool10 = bsxfun(@plus,Percentiles.Cool10,Perc_tmp4);
Percentiles.Cool10 = bsxfun(@plus,Percentiles.Cool10,Perc_tmp5);
Percentiles.Cool10 = bsxfun(@plus,Percentiles.Cool10,Perc_tmp6);
Percentiles.Cool10 = bsxfun(@ldivide,Percentiles.Amount,Percentiles.Cool10);
Percentiles.Cool10(Percentiles.Outside == -9999) = -9999;
clear Perc_tmp*

%25%
display('combined 25%')
Perc_tmp1 = Percentiles.Outside;
Perc_tmp2 = Percentiles.Outside;
Perc_tmp3 = Percentiles.Outside;
Perc_tmp4 = Percentiles.Outside;
Perc_tmp5 = Percentiles.Outside;
Perc_tmp6 = Percentiles.Outside;
Perc_tmp1(Percentiles.Model_1 == 2 |Percentiles.Model_1 == 1) = 1; 
Perc_tmp2(Percentiles.CN == 2 | Percentiles.CN == 1) = 1; 
Perc_tmp3(Percentiles.LPJ == 2 |Percentiles.LPJ == 1) = 1; 
Perc_tmp4(Percentiles.Benefit == 2 | Percentiles.Benefit == 1) = 1; 
Perc_tmp5(Percentiles.ScholesInt == 2 | Percentiles.ScholesInt == 1) = 1; 
Perc_tmp6(Percentiles.ScholesSA == 2 | Percentiles.ScholesSA == 1 ) = 1; 
Percentiles.Cool25 = bsxfun(@plus,Perc_tmp1,Perc_tmp2);
Percentiles.Cool25 = bsxfun(@plus,Percentiles.Cool25,Perc_tmp3);
Percentiles.Cool25 = bsxfun(@plus,Percentiles.Cool25,Perc_tmp4);
Percentiles.Cool25 = bsxfun(@plus,Percentiles.Cool25,Perc_tmp5);
Percentiles.Cool25 = bsxfun(@plus,Percentiles.Cool25,Perc_tmp6);
Percentiles.Cool25 = bsxfun(@ldivide,Percentiles.Amount,Percentiles.Cool25);
Percentiles.Cool25(Percentiles.Outside == -9999) = -9999;
clear Perc_tmp*

%90%
display('combined 90%')
Perc_tmp1 = Percentiles.Outside;
Perc_tmp2 = Percentiles.Outside;
Perc_tmp3 = Percentiles.Outside;
Perc_tmp4 = Percentiles.Outside;
Perc_tmp5 = Percentiles.Outside;
Perc_tmp6 = Percentiles.Outside;
Perc_tmp1(Percentiles.Model_1 == 5) = 1; 
Perc_tmp2(Percentiles.CN == 5) = 1; 
Perc_tmp3(Percentiles.LPJ == 5) = 1; 
Perc_tmp4(Percentiles.Benefit == 5) = 1; 
Perc_tmp5(Percentiles.ScholesInt == 5) = 1; 
Perc_tmp6(Percentiles.ScholesSA == 5) = 1; 
Percentiles.Warm90 = bsxfun(@plus,Perc_tmp1,Perc_tmp2);
Percentiles.Warm90 = bsxfun(@plus,Percentiles.Warm90,Perc_tmp3);
Percentiles.Warm90 = bsxfun(@plus,Percentiles.Warm90,Perc_tmp4);
Percentiles.Warm90 = bsxfun(@plus,Percentiles.Warm90,Perc_tmp5);
Percentiles.Warm90 = bsxfun(@plus,Percentiles.Warm90,Perc_tmp6);
Percentiles.Warm90 = bsxfun(@ldivide,Percentiles.Amount,Percentiles.Warm90);
Percentiles.Warm90(Percentiles.Outside == -9999) = -9999;
clear Perc_tmp*

%75%
display('combined 75%')
Perc_tmp1 = Percentiles.Outside;
Perc_tmp2 = Percentiles.Outside;
Perc_tmp3 = Percentiles.Outside;
Perc_tmp4 = Percentiles.Outside;
Perc_tmp5 = Percentiles.Outside;
Perc_tmp6 = Percentiles.Outside;
Perc_tmp1(Percentiles.Model_1 == 4 |Percentiles.Model_1 == 5) = 1; 
Perc_tmp2(Percentiles.CN == 4 | Percentiles.CN == 5) = 1; 
Perc_tmp3(Percentiles.LPJ == 4 |Percentiles.LPJ == 5) = 1; 
Perc_tmp4(Percentiles.Benefit == 4 | Percentiles.Benefit == 5) = 1; 
Perc_tmp5(Percentiles.ScholesInt == 4 | Percentiles.ScholesInt == 5) = 1; 
Perc_tmp6(Percentiles.ScholesSA == 4 | Percentiles.ScholesSA == 5) = 1; 
Percentiles.Warm75 = bsxfun(@plus,Perc_tmp1,Perc_tmp2);
Percentiles.Warm75 = bsxfun(@plus,Percentiles.Warm75,Perc_tmp3);
Percentiles.Warm75 = bsxfun(@plus,Percentiles.Warm75,Perc_tmp4);
Percentiles.Warm75 = bsxfun(@plus,Percentiles.Warm75,Perc_tmp5);
Percentiles.Warm75 = bsxfun(@plus,Percentiles.Warm75,Perc_tmp6);
Percentiles.Warm75 = bsxfun(@ldivide,Percentiles.Amount,Percentiles.Warm75);
Percentiles.Warm75(Percentiles.Outside == -9999) = -9999;
clear Perc_tmp*

%In between%
display('In between')
Perc_tmp1 = Percentiles.Outside;
Perc_tmp2 = Percentiles.Outside;
Perc_tmp3 = Percentiles.Outside;
Perc_tmp4 = Percentiles.Outside;
Perc_tmp5 = Percentiles.Outside;
Perc_tmp6 = Percentiles.Outside;
Perc_tmp1(Percentiles.Model_1 == 3) = 1; 
Perc_tmp2(Percentiles.CN == 3) = 1; 
Perc_tmp3(Percentiles.LPJ == 3) = 1; 
Perc_tmp4(Percentiles.Benefit == 3) = 1; 
Perc_tmp5(Percentiles.ScholesInt == 3) = 1; 
Perc_tmp6(Percentiles.ScholesSA == 3) = 1; 
Percentiles.InBetween = bsxfun(@plus,Perc_tmp1,Perc_tmp2);
Percentiles.InBetween = bsxfun(@plus,Percentiles.InBetween,Perc_tmp3);
Percentiles.InBetween = bsxfun(@plus,Percentiles.InBetween,Perc_tmp4);
Percentiles.InBetween = bsxfun(@plus,Percentiles.InBetween,Perc_tmp5);
Percentiles.InBetween = bsxfun(@plus,Percentiles.InBetween,Perc_tmp6);
Percentiles.InBetween = bsxfun(@ldivide,Percentiles.Amount,Percentiles.InBetween);
Percentiles.InBetween(Percentiles.Outside == -9999) = -9999;
clear Perc_tmp*
%% save data as ascii
Cool10 = Percentiles.Cool10;
save('Cool10_grazing.asc', 'Cool10', '-ascii')
clear Cool10
Cool25 = Percentiles.Cool25;
save('Cool25_grazing.asc', 'Cool25', '-ascii')
clear Cool25
Warm75 = Percentiles.Warm75;
save('Warm75_grazing.asc', 'Warm75', '-ascii')
clear Warm75
Warm90 = Percentiles.Warm90;
save('Warm90_grazing.asc', 'Warm90', '-ascii')
clear Warm90
InBetween = Percentiles.InBetween;
save('InBetween_grazing.asc', 'InBetween', '-ascii')
clear InBetween