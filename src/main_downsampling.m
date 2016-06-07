% Format precision
format long;

% ######## Parameters ########
nb_trials_training = 27; %11; % Ratio train/test for getting an accuracy score
highpass_filter = 1; %Hz
notch_filter = [49 51]; %Hz
path = 'H:\MasterProject\RawEEGData\actiCAP\SAL01\Online\';
name = 'motor-imagery-csp-4-online-[2015.11.01-11.54.07].vhdr';
hand = 1;
k = 3;
level = 5;
wavelet = 'sym1';
classifier = 1;
epoch_interval = [0 1];
% ############################

%Random seed
rand('state',sum(100*clock));

count = 1;

% Downsampling with different sampling rate
for i = (2000:-100:200)
    display(['#### Sampling rate ', num2str(i), ' - eeg variable n°', ...
        num2str(count)]);
    
    % #### Load the data file into EEGLAB toolbox
    EEG = pop_loadbv(path, name);
    EEG.setname = 'raw';
        
    alleeg_ds(count) =  pop_resample(EEG, i);
    count = count + 1;
end


for i = (1: length(alleeg_ds))
    display(['#### dataset ', num2str(i)]);
    start = tic();
    
    % #### Preprocessing
    alleeg_ds2(1) = alleeg_ds(i);
    
    % Re-reference (average) ---------> 23.8934%
    EEG = pop_reref( alleeg_ds(i), []);
    EEG.setname = 'Alistair re-referencing average';
    alleeg_ds2(2) = EEG;
    
    % 2.2 High pass filter (1Hz) ---------> 20.9124%
    %EEG = pop_eegfiltnew(EEG, [], highpass, 424, true, [], 0);
    EEG = pop_eegfiltnew(EEG, [], highpass_filter, 6600, true, [], 0);

    EEG.setname = 'Alistair highpass filter';
    alleeg_ds2(3) = EEG;
    
    % 2.3 Notch filter ------------> 20.8672%
    %EEG = pop_eegfiltnew(EEG, notch(1), notch(2), 212, 1, [], 0);
    EEG = pop_eegfiltnew(EEG, notch_filter(1), notch_filter(2), 6600, 1, [], 0);
    EEG.setname = 'Alistair notch filter';
    alleeg_ds2(4) = EEG;
    
    EEG = pop_autobsseog( EEG, [128], [128], 'sobi', {'eigratio', [1000000]}, 'eog_fd', {'range',[1  5]});
    alleeg_ds2(5) = EEG;
    
    % #### X: Epoch extraction
    EEG = pop_epoch( EEG, {'S769' 'S770'}, epoch_interval, 'newname', ...
        'Alistair epochs', 'epochinfo', 'yes');
    alleeg_ds2(6) = EEG;

    % Id of the latest preprocessed eeglab variable
    nb_dataset = length(alleeg_ds2);

    % #### 1: convert the data loaded on eeglab into a csv file and mat variables
    convertEEGLABdataIntoCsv(nb_dataset, alleeg_ds2, hand);
    load('data_events.mat');

    % #### 2: Features extraction
    mat_features = featuresExtraction(final_mat_X, level, ..., 
        alleeg_ds2(nb_dataset).trials, wavelet, alleeg_ds2(nb_dataset).nbchan);

    % #### 3: Features selection
    [bestMat, ~, ann_net] = featureSelection( mat_features, 20, size(mat_features, 2), ...
            200, k, 40, nb_trials_training, ex_events, 0.8, 0.1, classifier);

    % #### X: Apply k-NN
    disp(['###################### Final accuracy ############################']);
    [predictions, accuracy, ~] = classifiers(classifier, bestMat, ex_events,...
        nb_trials_training, alleeg_ds2(nb_dataset).trials, k, ann_net);
    
    times(i) = toc(start);
    accuracies(i) = accuracy;
end

times
accuracies