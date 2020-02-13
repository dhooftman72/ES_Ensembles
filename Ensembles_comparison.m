% Comparison function between Rho's and Deviances of Ensebles and Models
% Acuuracy vs Uncertainty regressions
clear all
clc
load('Ensembles.mat') % the accuracies of models (from Willcock et al. 2019) and Ensembles & Inversed deviances per point (calculated in Makeresults.m)
%% Preset admin
Comparisons = dataset({'Dummy'},'Varnames',char('Comparison'));
Original_Means = dataset({'Dummy'},'Varnames',char('Comparison'));
Best_model = dataset({'Dummy'},'Varnames','Service');
EnsembleArray = {'Model','Model','Mean','Median','Complexity','Deviance_Scaled','Rho_Scaled','Model'};
EnsembleArrayTxt = {'Best_Model','Models','Mean','Median','Complexity','Deviance','Rho','ModelAllInd'};
ImprovementValues.Default.Yrho = dataset(NaN,'Varnames',(genvarname({[char(EnsembleArrayTxt(1)),'vs',char(EnsembleArrayTxt(2))]})));
ImprovementValues.Full.Yrho = dataset(NaN,'Varnames',(genvarname({[char(EnsembleArrayTxt(length(EnsembleArray))),'vs',char(EnsembleArrayTxt(8))]})));
ImprovementValues.Default.Ydevi = dataset(NaN,'Varnames',(genvarname({[char(EnsembleArrayTxt(1)),'vs',char(EnsembleArrayTxt(2))]})));
ImprovementValues.Full.Ydevi = dataset(NaN,'Varnames',(genvarname({[char(EnsembleArrayTxt(length(EnsembleArray))),'vs',char(EnsembleArrayTxt(8))]})));
comparison = 1; % set array label (elevated in line 188)
%% All comparisons
for i = 1:length(EnsembleArray)
    % Pick the combinations
    jes = (i+1):(length(EnsembleArray)-1);
    AllInd = 0;
    if i == 1
        jes = (i+1):(length(EnsembleArray));
    elseif i == length(EnsembleArray)
        AllInd = 1;
        jes = 3:(length(EnsembleArray)-1);
    end
    for j = jes
        if j == 8 && i == 1 % only when i == 1
           AllInd = 1;
        end
        if strcmp('Model',EnsembleArray(i)) == 1 
            A = find(strcmp(EnsembleArray(i),Ensembles.ModelType)==1);
        else
            A =  find(strcmp('Mean',Ensembles.EnsembleType)==1);
        end
        if strcmp('Model',EnsembleArray(j)) == 1
            B = find(strcmp(EnsembleArray(j),Ensembles.ModelType)==1);
        else
            B = find(strcmp(EnsembleArray(j),Ensembles.EnsembleType)==1);
        end
        Comparisons.Comparison(comparison,1) = {[char(EnsembleArrayTxt(i)),'vs',char(EnsembleArrayTxt(j))]};
        Comparisons.A(comparison,1) = {char(EnsembleArrayTxt(i))};
        Comparisons.B(comparison,1) = {char(EnsembleArrayTxt(j))};
        Original_Means.Comparison = Comparisons.Comparison;
        Original_Means.A = Comparisons.A;
        Original_Means.B = Comparisons.B;
        %% the actual comparisons
        if strcmp('Best_Model',EnsembleArrayTxt(i)) ~= 1
            Y_A = cell2mat(Ensembles.Rho([A]));
            Ydevi_A = cell2mat(Ensembles.Deviation([A]));
            Y_B = cell2mat(Ensembles.Rho([B]));
            Ydevi_B = cell2mat(Ensembles.Deviation([B]));
            Sets = (Ensembles.UniqueSetID([A]));
            Sets_B = (Ensembles.UniqueSetID([B]));
            mean_YA = (Y_A);
            mean_YB = (Y_B);
            mean_YAdevi =(Ydevi_A);
            mean_YBdevi= (Ydevi_B);
            counts = 1;
            for t = 1:max(Sets)
                if AllInd == 0; % with averaging method
                    Y_A2 =  ((Y_A(Sets == t)+1)./2);
                    Y_B2 = nanmean(((Y_B(Sets_B == t)+1)./2));
                    Ycor(t) = nanmean(Y_B2./Y_A2); %#ok<*SAGROW
                    clear Y_A2 Y_B2
                    Y_A2 = Ydevi_A(Sets == t);
                    Y_B2 = nanmean(Ydevi_B(Sets_B == t));
                    Ydevicor(t) = nanmean(Y_B2./Y_A2); %#ok<*SAGROW>
                    clear Y_A2 Y_B2
                else % with full individual comparisons
                    Y_A2 =  ((Y_A(Sets == t)+1)./2);
                    counttmp = counts + length(Y_A2);
                    Y_B2 = nanmean(((Y_B(Sets_B == t)+1)./2));
                    Ycor(counts:(counttmp-1),1) = (Y_B2./Y_A2);
                    clear Y_A2 Y_B2
                    Y_A2 = Ydevi_A(Sets == t);
                    Y_B2 = nanmean(Ydevi_B(Sets_B == t));
                    Ydevicor(counts:(counttmp-1),1) = (Y_B2./Y_A2);
                    clear Y_A2 Y_B2
                    counts = counttmp;
                end
            end
        else  % Best models
            Y_A = cell2mat(Ensembles.Rho([A]));
            Ydevi_A = cell2mat(Ensembles.Deviation([A]));
            Model_s = (Ensembles.Model([A]));
            Y_B = cell2mat(Ensembles.Rho([B]));
            Ydevi_B = cell2mat(Ensembles.Deviation([B]));
            Sets = (Ensembles.UniqueSetID([A]));
            Sets_B = (Ensembles.UniqueSetID([B]));
            mean_YB = (Y_B);
            mean_YBdevi= (Ydevi_B);
            counts = 1;
            for t = 1:max(Sets)
                if  AllInd == 0; % with averaging method
                    Modelset = Model_s(Sets == t);
                    Y_A2 =  ((Y_A(Sets == t)+1)./2);
                    Y_B2 = nanmean(((Y_B(Sets_B == t)+1)./2),1);
                    Ycor(t) = max(Y_A2)./Y_B2;
                    mean_YA(t) = max(Y_A(Sets == t));
                    testMax = find(Y_A2 == max(Y_A2));
                    if strcmp('Mean',EnsembleArrayTxt(j)) ~= 1
                        if length(testMax) > 1
                            nameBest = ['joined_',char(Modelset(testMax(1))),'_',char(Modelset(testMax(2)))];
                        else
                            nameBest = char(Modelset(testMax));
                        end
                        Best_model.Service(t,1) = Ensembles.Service([B(t)]);
                        Best_model.Validation(t,1) = Ensembles.Validation([B(t)]);
                        Best_model.Rho(t,1) = {nameBest};
                    end
                    clear Y_A2 Y_B2 ycons testMax name
                    
                    Y_A2 =  Ydevi_A(Sets == t);
                    Y_B2 = nanmean(Ydevi_B(Sets_B == t),1);
                    Ydevicor(t) = max(Y_A2)./Y_B2;
                    mean_YAdevi(t) = max(Y_A2);
                    testMax = find(Y_A2 == max(Y_A2));
                    if strcmp('Mean',EnsembleArrayTxt(j)) ~= 1
                        if length(testMax) > 1
                            nameBest = ['joined_',char(Modelset(testMax(1))),'_',char(Modelset(testMax(2)))];
                        else
                            nameBest = char(Modelset(testMax));
                        end
                        Best_model.Devi(t,1) = {nameBest};
                    end
                    clear Y_A2 Y_B2 name
                else % with full individual comparisons
                    Y_A2 =  ((Y_A(Sets == t)+1)./2);
                    counttmp = counts + length(Y_A2);
                    Modelset = Model_s(Sets == t);
                    Y_B2 = ((Y_B(Sets_B == t)+1)./2);
                    Ycor(counts:(counttmp-1),1) = max(Y_A2)./Y_B2;
                    clear Y_A2 Y_B2
                    mean_YA(t) = max(Y_A(Sets == t));
                    Y_A2 =  Ydevi_A(Sets == t);
                    Y_B2 = (Ydevi_B(Sets_B == t));
                    Ydevicor(counts:(counttmp-1),1) = max(Y_A2)./Y_B2; %#ok<*SAGROW>
                    mean_YAdevi(t) = max(Y_A2);
                    clear Y_A2 Y_B2 name
                    counts = counttmp;
                end
            end
        end
        Y(:,1) = Ycor;
        Ydevi(:,1) = Ydevicor;
        clear t
        clear Ycor Ydevicor
        if AllInd == 0
            ImprovementValues.Default.Yrho.(genvarname(char(Comparisons.Comparison(comparison,1)))) = Y;
            ImprovementValues.Default.Ydevi.(genvarname(char(Comparisons.Comparison(comparison,1)))) = Ydevi;
        else
            ImprovementValues.Full.Yrho.(genvarname(char(Comparisons.Comparison(comparison,1)))) = Y;
            ImprovementValues.Full.Ydevi.(genvarname(char(Comparisons.Comparison(comparison,1)))) = Ydevi;
        end
        [~,Pnormal] = kstest(Y); %test for normality
        [~,Pvalue,~,stats] = ttest(Y,1);
        Comparisons.Significance(comparison,1) = Pvalue;
        Comparisons.DF(comparison,1) = stats.df;
        Comparisons.Tstat(comparison,1) = stats.tstat;
        Comparisons.PNormal(comparison,1) = Pnormal;
        Comparisons.meanB(comparison,1) = ((nanmean(Y(1:length(Y))))-1);
        Comparisons.SEM(comparison,1) = (nanstd(Y(1:length(Y)))./sqrt(16));
        clear Pvalue Y stats

        [~,Pnormal] = kstest(Ydevi); %test fort normality
        [~,Pvalue,~,stats] = ttest(Ydevi,1);
        Comparisons.Significance_Devi(comparison,1) = Pvalue;
        Comparisons.DF_Devi(comparison,1) = stats.df;
        Comparisons.Tstat_Devi(comparison,1) = stats.tstat;
        Comparisons.PNormal_Devi(comparison,1) = Pnormal;
        Comparisons.meanB_Devi(comparison,1) = ((nanmean(Ydevi(1:length(Ydevi))))-1);
        Comparisons.SEM_Devi(comparison,1) = (nanstd(Ydevi(1:length(Ydevi)))./sqrt(16));
        clear Pvalue Ydevi stats
        
        Original_Means.MeanAorg(comparison,1) = nanmean(mean_YA);
        Original_Means.MeanBorg(comparison,1) = nanmean(mean_YB);
        Original_Means.SemAorg(comparison,1) = (nanstd(mean_YA))./sqrt(16);
        Original_Means.SemBorg(comparison,1) = (nanstd(mean_YB))./sqrt(16);
        Original_Means.MeanAorgDevi(comparison,1) =  nanmean(mean_YAdevi);
        Original_Means.MeanBorgDevi(comparison,1) = nanmean(mean_YBdevi);
        Original_Means.SemAorgDevi(comparison,1) =  (nanstd(mean_YAdevi))./sqrt(16);
        Original_Means.SemBorgDevi(comparison,1) =  (nanstd(mean_YBdevi))./sqrt(16);
        
        Comparisons.DifSEM(comparison,1) = (Comparisons.SEM(comparison,1)/Comparisons.meanB(comparison,1)).*...
            abs((Original_Means.MeanAorg(comparison,1)-Original_Means.MeanBorg(comparison,1)));
        Comparisons.DifSEMdevi(comparison,1) = (Comparisons.SEM_Devi(comparison,1)/Comparisons.meanB_Devi(comparison,1)).*...
            abs((Original_Means.MeanAorgDevi(comparison,1)-Original_Means.MeanBorgDevi(comparison,1)));
        
        comparison = comparison +1;
    end
end
save('Outputs','Comparisons','Deviations','ImprovementValues','Best_model', 'Original_Means');
clear all
load('Outputs.mat')
delete('Outputs.mat')

%% % Regressions of Accurcay vs Uncertainty
Service_nots = unique(Deviations.Service,'stable');
whichones = (1:length(Deviations.Service));%
Points = dataset(Deviations.Service(whichones),'Varnames','Services');
Points.SEM = Deviations.MeanStd (whichones);
CorFac = Deviations.N(whichones);
Points.SEM = Points.SEM./sqrt(CorFac);
Points.Means = Deviations.Mean(whichones);
Points.MeanLogit = log10((Points.Means./(1-Points.Means))+1);
Points.MeanLogit(Points.MeanLogit>0.99) = 0.99;
% Interaction model
[~,outs,stats] = anovan(Points.MeanLogit,{Points.Services,Points.SEM},'sstype',2,...
    'model','full','continuous', 2,'display', 'off',...
    'varnames', {'Service','Std'});
Points.MeanLogitCorrected =  Points.MeanLogit- stats.resid;
Regressions.Mean = outs;
Regressions.EffectPerES = dataset(stats.coeffs((length(stats.coeffs))-(5)),'Varnames',char(Service_nots(1)));
for i = 5:-1:1
    Regressions.EffectPerES.(genvarname(char(Service_nots(7-i)))) = stats.coeffs((length(stats.coeffs))-(i-1));
end
% Ome-way model for R2 only
X = Points.SEM;
Y = Points.MeanLogit;
ds = dataset(Y,X);
ds.Services = ordinal(Points.Services);
mdl = LinearModel.fit(ds,' Y ~ X + Services');
anova(mdl,'component',1);
Regressions.Mean(1,8) = {'R2'};
Regressions.Mean(2,8) = {mdl.Rsquared.Adjusted};
clear outs stats Y X
%%
save('Outputs','Comparisons','ImprovementValues','Best_model', 'Original_Means','Regressions','Points');
clear all
load('Outputs.mat')
