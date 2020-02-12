% Comparison function between Rho's and Deviances of Ensebles and Models
% Acuuracy vs Uncertainty regressions
clear all
clc
load('Ensembles.mat')
warning off
Comparisons = dataset({'Dummy'},'Varnames',char('Comparison'));
Original_Means = dataset({'Dummy'},'Varnames',char('Comparison'));
Best_model = dataset({'Dummy'},'Varnames','Service');
ImprovementValues.Default.Yrho = dataset(NaN,'Varnames',(genvarname(char({'ModelvsMean'}))));
ImprovementValues.Full.Yrho = dataset(NaN,'Varnames',(genvarname(char({'ModelvsMean'}))));
ImprovementValues.Default.Ydevi = dataset(NaN,'Varnames',(genvarname(char({'ModelvsMean'}))));
ImprovementValues.Full.Ydevi = dataset(NaN,'Varnames',(genvarname(char({'ModelvsMean'}))));
for comparison = 1:27
    if comparison == 1 || comparison == 22
        A = find(strcmp('Model',Ensembles.ModelType)==1);
        B = find(strcmp('Mean',Ensembles.EnsembleType)==1);
        Comparisons.Comparison(comparison,1) = {'ModelvsMean'} ;
        Comparisons.A(comparison,1) = {'Models'};
        Comparisons.B(comparison,1) = {'Mean'};
    elseif comparison ==2 || comparison == 23
        A = find(strcmp('Model',Ensembles.ModelType)==1);
        B = find(strcmp('Median',Ensembles.EnsembleType)==1);
        Comparisons.Comparison(comparison,1) = {'ModelvsMedian'} ;
        Comparisons.A(comparison,1) = {'Models'};
        Comparisons.B(comparison,1) = {'Median'};
    elseif comparison == 3 || comparison == 24
        A = find(strcmp('Model',Ensembles.ModelType)==1);
        B = find(strcmp('Complexity',Ensembles.EnsembleType)==1);
        Comparisons.Comparison(comparison,1) = {'ModelvsComplexity'} ;
        Comparisons.A(comparison,1) = {'Models'};
        Comparisons.B(comparison,1) = {'Complexity'};
    elseif comparison == 4 || comparison == 25
        A = find(strcmp('Model',Ensembles.ModelType)==1);
        B = find(strcmp('Deviance_Scaled',Ensembles.EnsembleType)==1);
        Comparisons.Comparison(comparison,1) = {'ModelvsScaledDeviance'} ;
        Comparisons.A(comparison,1) = {'Models'};
        Comparisons.B(comparison,1) = {'Deviance'};
    elseif comparison == 5 || comparison == 26
        A = find(strcmp('Model',Ensembles.ModelType)==1);
        B = find(strcmp('Rho_Scaled',Ensembles.EnsembleType)==1);
        Comparisons.Comparison(comparison,1) = {'ModelvsScaledRho'} ;
        Comparisons.A(comparison,1) = {'Models'};
        Comparisons.B(comparison,1) = {'Rho'};
        
        % One by one comparisons
    elseif comparison == 6 %
        A =  find(strcmp('Mean',Ensembles.EnsembleType)==1);
        B = find(strcmp('Median',Ensembles.EnsembleType)==1);
        Comparisons.Comparison(comparison,1) = {'MeanvsMedian'};
        Comparisons.A(comparison,1) = {'Mean'};
        Comparisons.B(comparison,1) = {'Median'};
    elseif comparison == 7
        A =  find(strcmp('Mean',Ensembles.EnsembleType)==1);
        B = find(strcmp('Complexity',Ensembles.EnsembleType)==1);
        Comparisons.Comparison(comparison,1) = {'MeanvsComplexity'};
        Comparisons.A(comparison,1) = {'Mean'};
        Comparisons.B(comparison,1) = {'Complexity'};
    elseif comparison == 8 %32
        A =  find(strcmp('Mean',Ensembles.EnsembleType)==1);
        B = find(strcmp('Deviance_Scaled',Ensembles.EnsembleType)==1);
        Comparisons.Comparison(comparison,1) = {'MeanvsDeviance'};
        Comparisons.A(comparison,1) = {'Mean'};
        Comparisons.B(comparison,1) = {'Deviance'};
    elseif comparison == 9 %33
        A =  find(strcmp('Mean',Ensembles.EnsembleType)==1);
        B = find(strcmp('Rho_Scaled',Ensembles.EnsembleType)==1);
        Comparisons.Comparison(comparison,1) = {'MeanvsRho'};
        Comparisons.A(comparison,1) = {'Mean'};
        Comparisons.B(comparison,1) = {'Rho'};
        % Median
    elseif comparison == 10
        A =  find(strcmp('Median',Ensembles.EnsembleType)==1);
        B = find(strcmp('Complexity',Ensembles.EnsembleType)==1);
        Comparisons.Comparison(comparison,1) = {'MedianvsComplexity'};
        Comparisons.A(comparison,1) = {'Median'};
        Comparisons.B(comparison,1) = {'Complexity'};
    elseif comparison == 11 %37
        A =  find(strcmp('Median',Ensembles.EnsembleType)==1);
        B = find(strcmp('Deviance_Scaled',Ensembles.EnsembleType)==1);
        Comparisons.Comparison(comparison,1) = {'MedianvsDeviance'};
        Comparisons.A(comparison,1) = {'Median'};
        Comparisons.B(comparison,1) = {'Deviance'};
    elseif comparison == 12 %38
        A =  find(strcmp('Median',Ensembles.EnsembleType)==1);
        B = find(strcmp('Rho_Scaled',Ensembles.EnsembleType)==1);
        Comparisons.Comparison(comparison,1) = {'MedianvsRho'};
        Comparisons.A(comparison,1) = {'Median'};
        Comparisons.B(comparison,1) = {'Rho'};
        % Deviance
    elseif comparison == 13
        A =  find(strcmp('Complexity',Ensembles.EnsembleType)==1);
        B = find(strcmp('Deviance_Scaled',Ensembles.EnsembleType)==1);
        Comparisons.Comparison(comparison,1) = {'ComplexityvsDeviance'};
        Comparisons.A(comparison,1) = {'Complexity'};
        Comparisons.B(comparison,1) = {'Deviance'};
    elseif comparison == 14
        A =  find(strcmp('Complexity',Ensembles.EnsembleType)==1);
        B = find(strcmp('Rho_Scaled',Ensembles.EnsembleType)==1);
        Comparisons.Comparison(comparison,1) = {'ComplexityvsRho'};
        Comparisons.A(comparison,1) = {'Complexity'};
        Comparisons.B(comparison,1) = {'Rho'};
    elseif comparison == 15 %48
        A =  find(strcmp('Deviance_Scaled',Ensembles.EnsembleType)==1);
        B = find(strcmp('Rho_Scaled',Ensembles.EnsembleType)==1);
        Comparisons.Comparison(comparison,1) = {'Deviance_ScaledvsRho'};
        Comparisons.A(comparison,1) = {'Deviance_Scaled'};
        Comparisons.B(comparison,1) = {'Rho'};
        % Rho
        
        % % Best
    elseif comparison == 16 || comparison == 27
        A = find(strcmp('Model',Ensembles.ModelType)==1);
        B = find(strcmp('Model',Ensembles.ModelType)==1);
        Comparisons.Comparison(comparison,1) = {'Bestvsmodels'};
        Comparisons.A(comparison,1) = {'Models'};
        Comparisons.B(comparison,1) = {'Best'};
    elseif comparison== 17 %57
        A = find(strcmp('Model',Ensembles.ModelType)==1);
        B = find(strcmp('Mean',Ensembles.EnsembleType)==1);
        Comparisons.Comparison(comparison,1) = {'BestvsMean'} ;
        Comparisons.A(comparison,1) = {'Mean'};
        Comparisons.B(comparison,1) = {'Best'};
    elseif comparison == 18 %58
        A = find(strcmp('Model',Ensembles.ModelType)==1);
        B = find(strcmp('Median',Ensembles.EnsembleType)==1);
        Comparisons.Comparison(comparison,1) = {'BestvsMedian'} ;
        Comparisons.A(comparison,1) = {'Median'};
        Comparisons.B(comparison,1) = {'Best'};
    elseif comparison == 19
        A = find(strcmp('Model',Ensembles.ModelType)==1);
        B = find(strcmp('Complexity',Ensembles.EnsembleType)==1);
        Comparisons.Comparison(comparison,1) = {'BestvsComplexity'} ;
        Comparisons.A(comparison,1) = {'Complexity'};
        Comparisons.B(comparison,1) = {'Best'};
    elseif comparison == 20 %60
        A = find(strcmp('Model',Ensembles.ModelType)==1);
        B = find(strcmp('Deviance_Scaled',Ensembles.EnsembleType)==1);
        Comparisons.Comparison(comparison,1) = {'BestvsScaledDeviance'};
        Comparisons.A(comparison,1) = {'Scaled_Deviance'};
        Comparisons.B(comparison,1) = {'Best'};
    elseif comparison == 21 % 61
        A = find(strcmp('Model',Ensembles.ModelType)==1);
        B = find(strcmp('Rho_Scaled',Ensembles.EnsembleType)==1);
        Comparisons.Comparison(comparison,1) = {'BestvsScaledRho'};
        Comparisons.A(comparison,1) = {'Scaled_Rho'};
        Comparisons.B(comparison,1) = {'Best'};
    end
    Original_Means.Comparison = Comparisons.Comparison;
    Original_Means.A = Comparisons.A;
    Original_Means.B = Comparisons.B;
    
    if exist('A') == 1
        if comparison ~=[16:21,27]
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
                if comparison <=21
                    Y_A2 =  ((Y_A(Sets == t)+1)./2);
                    Y_B2 = nanmean(((Y_B(Sets_B == t)+1)./2));
                    Ycor(t) = nanmean(Y_B2./Y_A2); %#ok<*SAGROW
                    clear Y_A2 Y_B2
                    Y_A2 = Ydevi_A(Sets == t);
                    Y_B2 = nanmean(Ydevi_B(Sets_B == t));
                    Ydevicor(t) = nanmean(Y_B2./Y_A2); %#ok<*SAGROW>
                    clear Y_A2 Y_B2
                else
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
            
            % Best
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
                if comparison <=21
                    Modelset = Model_s(Sets == t);
                    Y_A2 =  ((Y_A(Sets == t)+1)./2);
                    Y_B2 = nanmean(((Y_B(Sets_B == t)+1)./2),1);
                    Ycor(t) = max(Y_A2)./Y_B2;
                    mean_YA(t) = max(Y_A(Sets == t));
                    testMax = find(Y_A2 == max(Y_A2));
                    if comparison == 17
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
                    if comparison == 17
                        if length(testMax) > 1
                            nameBest = ['joined_',char(Modelset(testMax(1))),'_',char(Modelset(testMax(2)))];
                        else
                            nameBest = char(Modelset(testMax));
                        end
                        Best_model.Devi(t,1) = {nameBest};
                    end
                    clear Y_A2 Y_B2 name
                else
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
        if comparison< 22
         ImprovementValues.Default.Yrho.(genvarname(char(Comparisons.Comparison(comparison,1)))) = Y;
        else
           ImprovementValues.Full.Yrho.(genvarname(char(Comparisons.Comparison(comparison,1)))) = Y;
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
        
         if comparison< 22
         ImprovementValues.Default.Ydevi.(genvarname(char(Comparisons.Comparison(comparison,1)))) = Ydevi;
        else
           ImprovementValues.Full.Ydevi.(genvarname(char(Comparisons.Comparison(comparison,1)))) = Ydevi;
         end
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
    end
end
save('Outputs','Comparisons','Ensembles','Deviations','ImprovementValues','Best_model', 'Original_Means');
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