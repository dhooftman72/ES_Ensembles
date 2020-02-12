% Statistics function for Willcock et al.
function [datapoints, RHO_all,PVAL_all, mean_double_deviation,deviation_point, xes, yes] = Accuracy_statistics(test,make_log) 
%Input: test
%Outputs: 
% datapoints; RHO_all; PVAL_all; mean_double_deviation; deviation_point,
% xes, yes

% clean data set
test(isinf(test)==1) = NaN;
orginal_order = 1:1:length(test(:,1));
orginal_order = orginal_order';
a=find((isnan(test(:,1))==1));
b=find((isnan(test(:,2))==1));
a = [a;b];
test(a,:) = [];  %#ok<*FNDSB>
orginal_order(a,:) = [];
missing_order = a;
datapoints = size(test,1);
clear a b

%% Make log if applicable, but not for Ensembles
testVar = test;
if make_log == 1
    testVar = log10(test+1);
end
%% Correlation stats: Spearman: Overall Accuracy
[RHO_all,PVAL_all] = corr(testVar(:,2),testVar(:,1),'type','Spearman'); 
%% Inverse deviance against a 1:1 line
% Winzorisation
x_range = testVar(:,2)./(prctile(testVar(:,2),95));
y_range = testVar(:,1)./(prctile(testVar(:,1),95));
x_range(x_range>1) = 1;
y_range(y_range>1) = 1;
clear deviation_point
deviation_point= abs(y_range-x_range); %% Accuracy per point   
mean_double_deviation = 1- ((sum(abs(y_range-x_range)))/datapoints); %%Accuracy overall
% %% recreate x and y ranges in orginal order, for individual model point
% storage to create ensembles
clear xes yes
xes(orginal_order) = x_range;
xes(missing_order) = NaN;
yes(orginal_order) = y_range;
yes(missing_order) = NaN;
clear x_range y_range
clear orginal_order missing_order
clear test
end