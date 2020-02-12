function [Results,Validation_points,Model_points,Weighting] = MakeResults(datapoints, RHO_all,...
    PVAL_all, mean_double_deviation,xes, yes,data_set,data_set_max,Names,Supply_combine,...
    People_combine,Supply,Results,Validation_points,Model_points,Weighting,output_file)
%% All calculations
if data_set == 1
Results.Models = dataset(NaN,'Varnames',char({'Data_set'}));
end
if data_set <= data_set_max 
    % Output for individual models as check on Willcock et al 2019
    Weighting.devi(data_set) = mean_double_deviation;% for use in weighting
    Weighting.Rho(data_set) = RHO_all; % for use in weighting
    Results.Models.Data_set(data_set,1) = data_set;
    Results.Models.Datapoints(data_set,1) = datapoints;
    Results.Models.RHO(data_set,1) = (round(RHO_all*10000))./10000;
    Results.Models.PVal(data_set,1) =  (round(double(PVAL_all)*10000))./10000;
    Results.Models.Inversed_deviance(data_set,1) =   (round(mean_double_deviation*10000))./10000;
    % Collect Individual data points
    Validation_points(data_set,:) = xes;
    Model_points(data_set,:) = yes;
end
if data_set ==  data_set_max
    display ('Summing results')
    [Results, Deviations] = model_ensembles(Results,Names,Supply_combine,...
    People_combine,Validation_points,Model_points,Weighting,Supply);
save(output_file,'Results','Deviations')
end
delete('output_files.mat')
delete('parameters.mat')
delete('validation_set.mat')
end

%%
function [Results, Deviations] = model_ensembles(Results,Names,Supply_combine,...
    People_combine,Validation_points,Model_points,Weighting,Supply)
% Ensemble module
% Here the bespoke ensembles are created, tested and stored
display ('Combining models')
make_log = 0;
Deviations = dataset({'Dummy'},'Varnames',char('Datapoint_name'));
%% Preparations
Namescor = Names;
% Set data_set
if Supply == 1
    vals = Supply_combine;
elseif Supply ==0
    vals = People_combine;
end
set_max = length(vals);
validation_store = Validation_points(vals(1),:);
val_tmp = [];
for i=1:1:set_max
    val_tmp = [val_tmp; Model_points(vals(i),:)]; %#ok<*AGROW>
end
model_values_store = val_tmp;
clear val_tmp i 

% Remove missing x-data
list = find(isnan(validation_store) == 1);
validation_store(list) = [];
model_values_store(:,list) = [];
Namescor(list) = [];
%% Values for weighting
Deviances = Weighting.devi(vals);
Rhos_overview =  ((( Weighting.Rho(vals))+1)./2);
Deviances_cor_max = max(Deviances);
Deviances_cor_range = Deviances_cor_max  - min(Deviances);
Rhos_cor_max = max(Rhos_overview);
Rhos_cor_range = Rhos_cor_max  -min(Rhos_overview);
for i = 1:length(vals)
    Deviances_cor(i) = Deviances_cor_max - ((Deviances_cor_max - Deviances(i)).*(1./Deviances_cor_range)); %#ok<*SAGROW>
    Rhos_cor(i) = Rhos_cor_max - ((Rhos_cor_max - Rhos_overview(i)).*(1./Rhos_cor_range));
end
Deviances_cor(Deviances_cor < 0.25) = 0.25;
Rhos_cor(Rhos_cor < 0.25) = 0.25;
clear Deviances_cor_max Deviances_cor_range Rhos_cor_max Rhos_cor_range i vals Deviances Rhos_overview

Ensemble_Names = {'Mean';'Median';'Complexity_weighted';'Deviance_weighted';'Rho_weighted'};
Result_model_txt = {'Datapoints';'RHO';'PVAL';'Inversed_deviance'};
Results.Ensemble = dataset([1;1;1;1],'ObsNames', Result_model_txt,'Varnames',char(Ensemble_Names(1)));

%% Ensemble stats
%% Mean
validation = validation_store';
model_values = (nanmean(model_values_store,1))'; %The recalculation method
test = [model_values,validation]; 

[datapoints, RHO_all,PVAL_all, mean_double_deviation,deviation_point, ~, ~] = Accuracy_statistics(test, make_log);
Results.Ensemble.(genvarname(char(Ensemble_Names(1)))) = [datapoints;RHO_all;double(PVAL_all);mean_double_deviation]; %output
clear mean_double_deviation RHO_all PVAL_all datapoints

% Uncertainty vs accuracy outputs ONLY DONE WITH MEAN
 Var_std =  nanstd(model_values_store); % Uncertainty of Ensemble
 Var_distance_combined = 1- deviation_point; % Accuracy of Ensemble Mean
for i = 1:1:length(model_values) % number of models per datapoint
   Modelpoints(i) = length(find(isnan(model_values_store(:,i))~=1));
end
clear deviation_point validation model_values i
% clean out NaN
Var_M(isinf(Var_distance_combined)==1) = NaN;
Var_std(isinf(Var_std)==1) = NaN;
a=find((isnan(Var_distance_combined)==1));
b=find((isnan(Var_std)==1));
a = unique([a;b]);
Var_distance_combined(a) = [];
Var_std(a) = [];
Modelpoints(a) = [];
Names2 = Namescor;
Names2(a) = [];
clear a b
%Store
Deviations.Datapoint_name = Names2;
Deviations.Mean = reshape(Var_distance_combined,(length(Var_distance_combined)),1);
Deviations.MeanStd =  reshape(Var_std,(length(Var_std)),1);
Deviations.N = reshape(Modelpoints,(length(Modelpoints)),1);
clear Var_Mean Var_std Names2 Modelpoints

%% Median stats
validation = validation_store';
model_values = (nanmedian(model_values_store,1))'; %The recalculation method
test = [model_values,validation]; 
[datapoints, RHO_all,PVAL_all, mean_double_deviation,~, ~, ~] = Accuracy_statistics(test, make_log);
Results.Ensemble.(genvarname(char(Ensemble_Names(2)))) = [datapoints;RHO_all;double(PVAL_all);mean_double_deviation]; %output
clear model_values validation
clear mean_double_deviation RHO_all PVAL_all datapoints
%% Weighted Mean Complexity
validation = validation_store';
model_values = model_values_store';
for i = 1:1:length(model_values)
    nomi = 0;
    denomi = 0;
    for j = 1:1: set_max
        if isnan(model_values(i,j))~= 1
            weight_value = (Weighting.complex(j));
            nomi = nomi + (model_values(i,j).* weight_value);
            denomi = denomi +  weight_value;
        end
    end
    Var_combined_corrected(i) = nomi/denomi;
    clear nomi
    % see https://stackoverflow.com/questions/30383270/how-do-i-calculate-the-standard-deviation-between-weighted-measurements
    % and https://en.wikipedia.org/wiki/Reduced_chi-squared_statistic
    nomi = 0;
    for j = 1:1: set_max
        if isnan(model_values(i,j))~= 1
              weight_value = (Weighting.complex(j));
            nomi = nomi + (((model_values(i,j)-Var_combined_corrected(i))^2).*weight_value);
        end
    end
    Var_std(i) =  nomi/denomi;
end
test = [Var_combined_corrected',validation]; %#ok<*NASGU>
[datapoints, RHO_all,PVAL_all, mean_double_deviation,~, ~, ~] = Accuracy_statistics(test, make_log);
Results.Ensemble.(genvarname(char(Ensemble_Names(3)))) = [datapoints;RHO_all;double(PVAL_all);mean_double_deviation]; %output
clear Var_combined_corrected validation model_values
clear mean_double_deviation RHO_all PVAL_all datapoints
%% Weighted Mean Deviance Corrected
validation = validation_store';
model_values = model_values_store';
for i = 1:1:length(model_values)
    nomi = 0;
    denomi = 0;
    for j = 1:1: set_max
        if isnan(model_values(i,j))~= 1
            weight_value = Deviances_cor(j);
            nomi = nomi + (model_values(i,j).* weight_value);
            denomi = denomi +  weight_value;
        end
    end
    Var_combined_corrected(i) = nomi/denomi;
    clear nomi
    
    % see https://stackoverflow.com/questions/30383270/how-do-i-calculate-the-standard-deviation-between-weighted-measurements
    % and https://en.wikipedia.org/wiki/Reduced_chi-squared_statistic
    nomi = 0;
    for j = 1:1: set_max
        if isnan(model_values(i,j))~= 1
            weight_value = Deviances_cor(j);           
            nomi = nomi + (((model_values(i,j)-Var_combined_corrected(i))^2).*weight_value);
        end
    end
    Var_std(i) =  nomi/denomi;
end
test = [Var_combined_corrected',validation];
[datapoints, RHO_all,PVAL_all, mean_double_deviation,~, ~, ~] = Accuracy_statistics(test, make_log);
Results.Ensemble.(genvarname(char(Ensemble_Names(4)))) = [datapoints;RHO_all;double(PVAL_all);mean_double_deviation]; %output
clear Var_combined_corrected validation model_values
clear mean_double_deviation RHO_all PVAL_all datapoints
%% Weighted Mean Rho Corrected
validation = validation_store';
model_values = model_values_store';
for i = 1:1:length(model_values)
    nomi = 0;
    denomi = 0;
    for j = 1:1: set_max
        if isnan(model_values(i,j))~= 1
            weight_value = Rhos_cor(j);
            nomi = nomi + (model_values(i,j).* weight_value);
            denomi = denomi +  weight_value;
        end
    end
    Var_combined_corrected(i) = nomi/denomi;
    clear nomi
    
    % see https://stackoverflow.com/questions/30383270/how-do-i-calculate-the-standard-deviation-between-weighted-measurements
    % and https://en.wikipedia.org/wiki/Reduced_chi-squared_statistic
    nomi = 0;
    for j = 1:1: set_max
        if isnan(model_values(i,j))~= 1
            nomi = nomi + (((model_values(i,j)-Var_combined_corrected(i))^2).*weight_value);
        end
    end
    Var_std(i) =  nomi/denomi;
end
test = [Var_combined_corrected',validation];

[datapoints, RHO_all,PVAL_all, mean_double_deviation,~, ~, ~] = Accuracy_statistics(test, make_log);
Results.Ensemble.(genvarname(char(Ensemble_Names(5)))) = [datapoints;RHO_all;double(PVAL_all);mean_double_deviation]; %output
clear Var_combined_corrected validation model_values

end