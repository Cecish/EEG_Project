function accuracy = main()
clear all;
% Format precision
format long;

% ######## Parameters ########
highpass_filter = 0.3; %Hz
low_pass_filter = 60; %Hz
notch_filter = [48 52]; %[58 62]; %Hz
outcome_bad_channels = 0; %1: The user wants as an extra to do the data analysis
% by considering only the worst channels (according to channels selection
% outcome), 0 otherwise
do_plot = 1; %1: plot graphs 0: plot nothing
% ############################
manipFuns = dataManipFunctions; 
featureSelec = featureSelection;

% Letting the user decide of some settings in a console menu
% [hand, k, path_data_file, name_data_file, level, wavelet, classifier] = menu();
[hand, k, path_data_file, name_data_file, level, wavelet, classifier, str_device] = ...
   readParameterFile();

%Random seed
rng(sum(100*clock));

% #### 0: Preprocesing
% EEGLAB (include loading, preprocessing and extraction the EEG data)
alleeg = eeglab_script(path_data_file, name_data_file, highpass_filter, ...
    low_pass_filter, notch_filter, str_device);

% Id of the latest preprocessed eeglab variable
nb_dataset = length(alleeg);

% #### 1: convert the data loaded on eeglab into a csv file and mat variables
[final_mat_X, ex_events_Y, ex_events] = convertEEGLABdataIntoCsv(nb_dataset,...
    alleeg, hand);
% load('data_events.mat'); Uncomment if you have saved the data extracted
% in the convertEEGLABdataIntoCsv function into Matlab variabes
    
% #### 2: Min-Max normalisation
width = size(final_mat_X, 2);
height = size(final_mat_X, 1);

minmax_X = manipFuns.MinMaxNorm(final_mat_X, height, width);

% #### 3: Best channels subset selection
%By default, all the channels are considered for the analysis
desired_nb_channels = alleeg(nb_dataset).nbchan;
reduced_minmax_X = minmax_X;
selected_channels = 1:alleeg(nb_dataset).nbchan;
struct_mat{1} = reduced_minmax_X;
% desired_nb_channels = width*0.25; %The best 25% of all channels
desired_nb_channels = width*0.5; %The best 50% of all channels
% desired_nb_channels = width*0.75; %The best 75% of all channels
%In case the user commented the 3 lines above (and wants to consider all the channels)
if ~isequal(desired_nb_channels, alleeg(nb_dataset).nbchan)
    display([num2str(desired_nb_channels), 'selected for the analysis'])
    [reduced_minmax_X, selected_channels, worst_minmax_X] = featureSelec.chanReduction(...
        minmax_X, ex_events_Y, desired_nb_channels);
        struct_mat{1} = reduced_minmax_X;

        if (outcome_bad_channels)
            struct_mat{2} = worst_minmax_X;
        end
end

for i = (1: size(struct_mat, 2))
    % #### 4: Features extraction
    [mat_features, ratio_features] = featuresExtraction(struct_mat{i}, level, ...
        alleeg(nb_dataset).trials, wavelet, alleeg(nb_dataset).srate, ...
        alleeg(nb_dataset).pnts, desired_nb_channels*(2-i) + (i-1)*...
        (alleeg(nb_dataset).nbchan - desired_nb_channels ));

    % Randomly shuffle mat_features
    ordering = randperm(length(ex_events));
    mat_features2 = mat_features(ordering, :);
    ex_events2 = ex_events(ordering);

    % #### 5: Features selection
    [bestMat, best_features] = featureSelec.featSelection( mat_features2,...
        size(mat_features2, 2), k, ex_events2, classifier, do_plot);
    % best_features = geneticAlgorithm(mat_features2, ex_events2, classifier, k);

    % #### 5 bis: Analysis of the best features selected
    if do_plot
        manipFuns.analysisSelectedFeatures(best_features, desired_nb_channels, ratio_features);
    end
    
    % #### 6: Apply classifier
    disp('######### Final accuracy with k-folds cross validation ##############');
    [accuracy, CVMdl] = classifiers(classifier, bestMat, ex_events2, k);

    disp('######### Channels selected ##############');
    selected_channels

    disp('###################### Confusion Matrix ############################');
    if isequal(classifier, 3) %If the classifier is the MLP
        C= zeros(2, 2);
        
        predictions = sim(CVMdl, bestMat');
        for z = (1:40)
           if isequal(ex_events2(z), 0) && isequal(ex_events2(z), predictions(z))
               C(1, 1) = C(1, 1) + 1;
           elseif isequal(ex_events2(z), 0) && ~isequal(ex_events2(z), predictions(z))
               C(1, 2) = C(1, 2) + 1;
           elseif isequal(ex_events2(z), 1) && isequal(ex_events2(z), predictions(z))
               C(2, 1) = C(2, 1) + 1;
           elseif isequal(ex_events2(z), 1) && ~isequal(ex_events2(z), predictions(z))
               C(2, 2) = C(2, 2) + 1;
           end
        end
    else
        predictions = CVMdl.kfoldPredict;
        % figure('Name', 'Confusion Matrix')
        [C,order] = confusionmat(ex_events2, predictions');
        C
        order
        % plotconfusion(ex_events2, predictions');
        sensitivity = C(1,1) / (C(1,1) + C(2,1));
        specificity = C(2,2) / (C(2,2)+C(1,2));
        display(['Sensitivity = ', num2str(sensitivity)]);
        display(['Specificity = ', num2str(specificity)]);
    end

    disp('###################### Proportions ############################');
    nb_zeros = sum(ex_events2 == 0);
    disp(['0: ', num2str(nb_zeros/length(ex_events2)*100), '%']);
    disp(['1: ', num2str((length(ex_events2)-nb_zeros)/length(ex_events2)*100), '%']);
    % When I was evaluating my classifiers' performance by dividing the dataset
    % into a training and a testing set (no cross validation method)
    % nb_zeros = sum(ex_events2(1:length(ex_events2)*0.7) == 0);
    % disp(['0: ', num2str(nb_zeros/length(ex_events2)*0.7*100), '%']);
    % disp(['1: ', num2str((length(ex_events2)*0.7-nb_zeros)/length(ex_events2)*0.7*100), '%']);
end
end