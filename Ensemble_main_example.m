% Example of a Steering main Function for Willcock et al. Ecosystem service model ensembles
clear all
clc
%% join validation data with comparator
for validation_set = 1:1:x
    save('validation_set', 'validation_set')
    clear all
    load('validation_set.mat')
    load('Water_data.mat')
    make_log = 1;
    first = 1;
    if validation_set == 1
        data_set_max = 14;
        output_file = 'Results_ES1_ValidationSet1.mat'
        Supply = 1; %1 Supply; 0 = people
    elseif validation_set == 2
        ......
    end
for data_set = 1:1:data_set_max
    if first ~=1
        save('parameters', 'validation_set','data_set',...
            'make_log','Results','Validation_points','Model_points','data_set_max',...
            'Supply_combine', 'People_combine', 'Weighting','Supply')
        save ('output_files.mat', 'output_file')
        clear all
        clear all
        load('parameters.mat')
        load('ES1_Data.mat')
        load ('output_files.mat')
    end
    first = 0;
    if validation_set ==1
        %ValdiationSet 1
        Names = Site_names;
        make_log = 1;
        Supply_combine = [3,4,5,6,7,8];
        People_combine = [9,10,11,12,13,14];
        Weighting.complex = [17,125,10,36,2,24];
        bio_comparator_1 = (ValdiationSet1_validation_volume./ValdiationSet1_size(:,1));
        people_comparator_1 = (ValdiationSet1_demand(:,1)./ValdiationSet1_pop_size(:,1));
        % independents
        if data_set == 1
            test = [(ValdiationSet1_pop_size(:,1)./ValdiationSet1_size(:,1)),bio_comparator_1]; %i1 
        elseif data_set == 2
            test = [ValdiationSet1_mean_dens(:,1), people_comparator_1]; % i2
            %Model 1
        elseif data_set == 3
            test = [(ValdiationSet1(:,1)./ValdiationSet1_size(:,1)),bio_comparator_1]; % Supply
            %Model 2
        elseif data_set == 4
            test = [(ValdiationSet1(:,2)./ValdiationSet1_size(:,1)),bio_comparator_1]; % Supply
            % Model 3
        elseif data_set == 5
            test = [(ValdiationSet1(:,3)./ValdiationSet1_size(:,1)),bio_comparator_1]; % Supply
            %Model 4
        elseif data_set == 6
            test = [(ValdiationSet1(:,4)./ValdiationSet1_size(:,1)),bio_comparator_1]; % Supply
            %Model 5
        elseif data_set == 7
            test = [(ValdiationSet1(:,5)./ValdiationSet1_size(:,1)),bio_comparator_1]; % Supply
            %Model 6
        elseif data_set == 8
            test = [(ValdiationSet1(:,6)./ValdiationSet1_size(:,1)),bio_comparator_1]; % Supply
            %Model 1
        elseif data_set == 9
            test = [(ValdiationSet1(:,1)./ValdiationSet1_pop_size(:,1)),people_comparator_1]; %People
            %Model 2
        elseif data_set == 10
            test = [(ValdiationSet1(:,2)./ValdiationSet1_pop_size(:,1)),people_comparator_1]; %People
            % Model 3
        elseif data_set == 11
            test = [(ValdiationSet1(:,3)./ValdiationSet1_pop_size(:,1)),people_comparator_1]; %People
            %Model 4
        elseif data_set == 12
            test = [(ValdiationSet1(:,4)./ValdiationSet1_pop_size(:,1)),people_comparator_1]; %People
            %Model 5
        elseif data_set == 13
            test = [(ValdiationSet1(:,5)./ValdiationSet1_pop_size(:,1)),people_comparator_1]; %People
            %Model 6
        elseif data_set == 14
            test = [(ValdiationSet1(:,6)./ValdiationSet1_pop_size(:,1)),people_comparator_1]; %People
        else
            display('wrong model assigment')
            break
        end
    elseif validation_set ==2
        ......................
    end
%% Calculate
[datapoints, RHO_all,PVAL_all, mean_double_deviation,~, xes, yes] = Accuracy_statistics(test, make_log);
if data_set == 1
    Results = [];
    Validation_points = [];
    Model_points = [];
end
[Results,Validation_points,Model_points,Weighting] = MakeResults(datapoints, RHO_all,...
    PVAL_all, mean_double_deviation,xes, yes,data_set,data_set_max,Names,Supply_combine,...
    People_combine,Supply,Results,Validation_points,Model_points,Weighting,output_file);
end
end
clear all
