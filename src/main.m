% Format precision
format long;

% ######## Parameters ########
nb_trials_training = 27; %11; % Ratio train/test for getting an accuracy score
highpass_filter = 1; %Hz
notch_filter = [49 51]; %Hz
% ############################

% Letting the user decide of some settings in a console menu
[ hand, k, path_data_file, name_data_file, level, wavelet, classifier, ...
    device ] = menu();

%Random seed
rand('state',sum(100*clock));


% EEGLAB (include loading, preprocessing and extraction the EEG data)
alleeg = eeglab_script(path_data_file, name_data_file, highpass_filter, ...
    notch_filter, device);

% Update. For my tests. Can be commented and specified directly in the
% parameters section
nb_dataset = length(alleeg);

% #### 1: convert the data loaded on eeglab into a csv file and mat variables
convertEEGLABdataIntoCsv(nb_dataset, alleeg, hand);
load('data_events.mat');

% #### 2: Features extraction
mat_features = featuresExtraction(final_mat_X, level, ..., 
    alleeg(nb_dataset).trials, wavelet, alleeg(nb_dataset).nbchan);

% #### 3: Features selection
bestMat = featureSelection( mat_features, 20, size(mat_features, 2), ...
        200, k, 40, nb_trials_training, ex_events, 0.8, 0.1, classifier);

% #### X: Apply k-NN
disp(['###################### Final accuracy ############################']);
[predictions, accuracy] = classifiers(classifier, bestMat, ex_events,...
    nb_trials_training, alleeg(nb_dataset).trials, k);

