% Author: C�cile Riquart
% Date 16/05/2016, 17:22

% Converts the usefull data from the ALLEEG eeglab variable into a csv file
% Params:
%   - nb_dataset: the ALLEEG variable is a structure of EEG vraibles.
%   nb_dataset is the i of the EEG from ALLEEG that should be considered
%   - alleeg: ALLEEG variable generated by eeglab after its use
%   - hand: 'r'/'l': glove put on the right/left hand
function alleegToCsv(nb_dataset, alleeg, hand)

    % If the index specified as argument -> error
    if ((nb_dataset > length(alleeg)) || (nb_dataset < 1))
       error('Wrong dataset number: %d', nb_dataset);
    end

    export_3Ddata = alleeg(nb_dataset).data;
        
    %Transform 3D matrix into 2D matrix (trials are marged)
    final_mat_X = threeD2twoD(export_3Ddata, alleeg(nb_dataset).trials);
    
    %Build target field
    [ex_events, ex_events_Y] = buildTargetField(alleeg(nb_dataset), hand, length(final_mat_X));
    
    %Save current dataset into a csv file
    dlmwrite('dataseXY1.csv', [final_mat_X ex_events_Y], 'delimiter', ';');
    %dlmwrite('enormetest.txt', [final_mat_X ex_events_Y], 'delimiter', ',');
    
    % Save useful information into matlab variables
    save('data_events', 'final_mat_X', 'ex_events_Y', 'ex_events');
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
function [ex_events, ex_events_Y] = buildTargetField(eeg, hand, length_mat)
    switch hand
        case 'r'
            ex_events = extract_events(eeg, 'S770');
        case 'l'
            ex_events = extract_events(eeg, 'S769');
        otherwise
            error('The glove is either on the right (r) or left (l) hand.\n You said: %d', hand);
    end
    
    %Producing a colunm with the target field associating to each record
    ex_events_Y = rightDimTarget(length_mat, ex_events);
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

% Producing a colunm with the target field associating to each record 
% Params:
%   - length_mat: height of the 2D matrix storing the recordings
%   - events: array of events (right or left) for each trial
% Return: target field (column)
function last_col = rightDimTarget(length_mat, events)

    last_col = [];

    for i = 1:length(events)
        m = repmat(events(i), length_mat/length(events), 1);
        
        last_col = vertcat(last_col, m);
    end
end
