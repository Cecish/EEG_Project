% Format precision
format long;

% ######## Parameters ########
%nb_dataset = 5; % index of the EEG variable to consider among the ALLEEG global variable
hand = 'r';      % hand which wears the robotic glove (either 'r'-right or 'l'-left)
nb_trials_training = 27; %11; % Ratio train/test for getting an accuracy score
k = 1;  % Number of nearest neighbours
path_data_file = 'H:\MasterProject\RawEEGData\EPOC\acquisition\ToUse-Alistair\';
name_data_file = 'motor-imagery-csp-1-acquisition-[2016.03.18-16.55.10]_alistair.vdhr';
name_dataset = 'raw Alistair';
highpass_filter = 1; %Hz
notch_filter = [49 51]; %Hz
extract_epochs = [0 0.640]; %[starting the event 769 or 770 - end of the event]ms
level = 5; %For the discrete wavelet wavelet transform
wavelet = 'sym1'; %For the discrete wavelet wavelet transform
% ############################

%Random seed
%reset(RandStream.getDefaultStream,sum(100*clock)); 
rand('state',sum(100*clock));

% EEGLAB (include loading, preprocessing and extraction the EEG data)
alleeg = eeglab_script(path_data_file, name_data_file, name_dataset, ...
    highpass_filter, notch_filter, extract_epochs);

% Update. For my tests. Can be commented and specified directly in the
% parameters section
nb_dataset = length(alleeg);

% #### 1: convert the data loaded on eeglab into a csv file and mat variables
convertEEGLABdataIntoCsv(nb_dataset, alleeg, hand);
load('data_events.mat');

% #### 2: Features extraction
mat_features = featuresExtraction(final_mat_X, level, ..., 
    alleeg(nb_dataset).trials, wavelet);

% #### 3: Features selection
bestMat = featureSelection( mat_features, 20, size(mat_features, 2), ...
        200, 0.001, k, 40, nb_trials_training, events, 0.8, 0.1);

% #### X: Apply k-NN
%truc = kNN(mat_features, ex_events_Y, nb_trials_training, alleeg(nb_dataset).trials, k);
%%%%truc = kNN(mat_features, ty, nb_trials_training, alleeg(nb_dataset).trials, k);
truc = kNN(bestMat, ex_events, nb_trials_training, alleeg(nb_dataset).trials, k);

