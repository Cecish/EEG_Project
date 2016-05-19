% Author: C�cile Riquart
% Date 16/05/2016, 17:22

% Converts the usefull data from the ALLEEG eeglab variable into a csv file
% Params:
%   - nb_dataset: the ALLEEG variable is a structure of EEG vraibles.
%   nb_dataset is the i of the EEG from ALLEEG that should be considered
%   - alleeg: ALLEEG variable generated by eeglab after its use
%   - hand: 'r'/'l': glove put on the right/left hand
% Update: MinMax normalisation and 1-NN for checking the accuracy score
function alleegToCsv(nb_dataset, alleeg, hand)

    format long;

    % If the index specified as argument -> error
    if ((nb_dataset > length(alleeg)) || (nb_dataset < 1))
       error('Wrong dataset number: %d', nb_dataset);
    end

    export_3Ddata = alleeg(nb_dataset).data;
        
    %Transform 3D matrix into 2D matrix (trials are marged)
    final_mat = threeD2twoD(export_3Ddata, alleeg(nb_dataset).trials);
    
    %Build target field
    ex_events = buildTargetField(alleeg(nb_dataset), hand);
    
    %Build final dataset (including target field)
    final_mat = concatenate(final_mat, ex_events);

    %Save current dataset into a csv file
    dlmwrite('dataset1.csv', final_mat, 'delimiter', ';');
    
    %reduced version of this dataset by taking a random 1 in 10 of the rows
    
    %MinMaxNormalisation
    width = length(final_mat(1, :));
    height = length(final_mat);
    %width_n = length(ndataset(1, :));
    %height_n = length(ndataset);
        
    dataset_minmax = MinMaxNorm(final_mat, height, width);
    %dataset_minmax = MinMaxNorm(ndataset, height_n, width_n);
    
    %Accuracy
    accuracy = OneNearestNeighbour(dataset_minmax, height, width);
    %accuracy = OneNearestNeighbour(dataset_minmax, height_n, width_n);
    disp(['Accuracy with 1-NN = ', num2str(accuracy)]);
    
    %save('final_mat', 'final_mat');
end

% Transform 3D matrix into 2D matrix (trials are marged)
% Params:
%   - export_data: 3D data
%   - nb_trials: Number of trials
% Return: corresponding 2D matrix
function final_mat = threeD2twoD(export_data, nb_trials)
    final_mat = [];
    
    for i=1:nb_trials
       mat_temp = (export_data(:, :, i))';
       final_mat = vertcat(final_mat, mat_temp);
    end
end

% Build the target field column
% Params: 
%   - eeg: 2D matrix storing the EEG records
%   - hand: glove worn on the right ('r') or left ('l') hand
% Return events associated to each records (in array)
function ex_events = buildTargetField(eeg, hand)
    switch hand
        case 'r'
            ex_events = extract_events(eeg, 'S770');
        case 'l'
            ex_events = extract_events(eeg, 'S769');
        otherwise
            error('The glove is either on the right (r) or left (l) hand.\n You said: %d', hand);
    end
end

% Auxiliary function for building the target field column
% Params: 
%   - alleeg: 2D matrix storing the recordings
%   - state: code for identifying the right or left event
% Return: the target field column
function ex_events = extract_events(alleeg, state)

    temp = alleeg.event;
    
    for i = 1: length(alleeg.event)
        if (isequal(temp(1, i).type, state))
            ex_events(i) = 1;
        else
            ex_events(i) = 0;
        end
    end
end

% Concatenating the target field column (transformed to fit with the 2D 
% matrix dimensions) to the 2D matrix
% Params:
%   - mat: 2D matrix storing the recordings
%   - events: array of events (right or left) for each trial
% Return: new 2D matrix with the target field column added at the end
function final_mat = concatenate(mat, events)

    last_col = [];
    length(mat)
    length(events)

    for i = 1:length(events)
        m = repmat(events(i), length(mat)/length(events), 1);
        
        last_col = vertcat(last_col, m);
    end
    
    for i = (1: 3280)
        if (i <= 2624)
            trainingset(i, :) = mat(i, :);
            y1(i) = last_col(i);
        else
            validationset(i - 2624, :) = mat(i, :);
            y2(i - 2624) = last_col(i);
        end
    end
    Mdl = fitcknn(trainingset, y1);
    rloss = resubLoss(Mdl)
    pred = predict(Mdl, validationset);
    
    count = 0;
    for i=(1: length(pred))
        if isequal(pred(i), y2(i))
            count = count +1;
        end
    end
    perc = count / length(pred);
    perc
    
    final_mat = [mat last_col];
end