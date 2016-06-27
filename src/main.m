clear all;
% Format precision
format long;

% ######## Parameters ########
nb_trials_training = 28; %11; % Ratio train/test for getting an accuracy score
highpass_filter = 1; %Hz
notch_filter = [58 62]; %Hz
% ############################

% Letting the user decide of some settings in a console menu
[hand, k, path_data_file, name_data_file, level, wavelet, classifier] = menu();

%Random seed
rng(sum(100*clock));

% #### 0: Preprocesing
% EEGLAB (include loading, preprocessing and extraction the EEG data)
alleeg = eeglab_script(path_data_file, name_data_file, highpass_filter, ...
    notch_filter);

% Id of the latest preprocessed eeglab variable
nb_dataset = length(alleeg);

% #### 1: convert the data loaded on eeglab into a csv file and mat variables
convertEEGLABdataIntoCsv(nb_dataset, alleeg, hand);
load('data_events.mat');

% #### 2: Min-Max normalisation
width = size(final_mat_X, 2);
height = size(final_mat_X, 1);
manipFuns = dataManipFunctions; 
    
minmax_X = manipFuns.MinMaxNorm(final_mat_X, height, width);

% #### 3: Features extraction
mat_features = featuresExtraction(minmax_X, level, ...
    alleeg(nb_dataset).trials, wavelet, alleeg(nb_dataset).nbchan);

% Randomly shuffle mat_features
ordering = randperm(length(ex_events));
mat_features2 = mat_features(ordering, :);
ex_events2 = ex_events(ordering);

% #### 4: Features selection
[bestMat, best_features] = featureSelection( mat_features2, size(mat_features2, 2), k,...
     ex_events2, classifier);
% geneticAlgorithm(mat_features2, ex_events2, classifier, k);
% best_feature = geneticAlgorithm(mat_features2, ex_events2, k, classifier);

% #### 5: Apply classifier
disp('######### Final accuracy with k-folds cross validation ##############');
[accuracy, CVMdl] = classifiers(classifier, bestMat, ex_events2, k);

disp('###################### Confusion Matrix ############################');
predictions = CVMdl.kfoldPredict;
% figure('Name', 'Confusion Matrix')
[C,order] = confusionmat(ex_events2, predictions');
C
order
% plotconfusion(ex_events2, predictions');

disp('###################### Proportions ############################');
nb_zeros = sum(ex_events2 == 0);
disp(['0: ', num2str(nb_zeros/length(ex_events2)*100), '%']);
disp(['1: ', num2str((length(ex_events2)-nb_zeros)/length(ex_events2)*100), '%']);
% When I was evaluating my classifiers' performance by dividing the dataset
% into a training and a testing set (no cross validation method)
% nb_zeros = sum(ex_events2(1:length(ex_events2)*0.7) == 0);
% disp(['0: ', num2str(nb_zeros/length(ex_events2)*0.7*100), '%']);
% disp(['1: ', num2str((length(ex_events2)*0.7-nb_zeros)/length(ex_events2)*0.7*100), '%']);


