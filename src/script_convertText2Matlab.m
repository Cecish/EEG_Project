% Author: Cécile Riquart
% Date 15/05/2016, 10:34

% Reads a text file containing EEG raw data and converts it into a Matlab
% format + create a variable for associating the outcome to each record if
% the input file corresponds to a training file (known outcome)
% Params:
% - raw EEG data file
% - Matlab file to be output
% - Binary value: is it an acquistion(1) or training(0) data file?
% - Records corresponding to the right hand
% - Records corresponding to the left hand
% - The experiments where undertaken with the right hand ('r') or left hand
% ('l')
% 5 variables are created in Matlab:
% - Channels names
% - Sampling rate (frequency)
% - Time
% - signals values (expressed in microVolts)
% - outcome (0: targeted hand closed, 1: Targeted hand opened, -1 : none)
function convertText2Matlab(inputTextFile, outputMatlabFile, training,...
    rightTrialsFile, leftTrialsFile, left_or_right)

    format long; %To be sure that no information is lost during the conversion

    temp = importdata(inputTextFile, ';');
    
    switch training
        case 0 %Testing, outcome will be predict with a previously trained classifier
            outcome = [];
        case 1 %Training, outcome is known
            %Getting the records associated to an event (right or left hand)
            right = importdata(rightTrialsFile, ';');
            left = importdata(leftTrialsFile, ';');
            %Building the outcome variable
            outcome = buildOutcome(temp.data(:, 1), right.data(:, 1), ...
                left.data(:, 1), left_or_right);
        otherwise
            error('The third argument is either true (1, the input file is an acquisition data file) or false (0)\nYour parameter: %d', training);
    end
            
    channels_name = temp.textdata(2:length(temp.textdata)-1); %channels names
    nb_channels = (length(temp.textdata)-2); %number of channels (exclude sampling rate and time columns)
    sampling_rate = temp.data(1, length(temp.textdata)); %sampling rate
    time = temp.data(:, 1); %time
    eeg_values = (temp.data(:, 2:length(temp.textdata)-1)).'; %eeg values
    
    % ---------------------- Variable as structure ------------------------
    %rawEEG(1).channels_name = temp.textdata(2:length(temp.textdata)-1); %channels names
    %rawEEG(1).nbchan = (length(temp.textdata)-2); %number of channels (exclude sampling rate and time columns)
    %rawEEG(1).srate = temp.data(1, length(temp.textdata)); %sampling rate
    %rawEEG(1).time = temp.data(:, 1); %time
    %rawEEG(1).data = temp.data(:, 2:length(temp.textdata)-1); %ee values
    %rawEEG(1).outcome = outcome;  %outcome
    % ---------------------------------------------------------------------
    
    %Save specified variables in outputMatlabFile
    save(outputMatlabFile, 'channels_name', 'nb_channels', 'sampling_rate', ...
        'time', 'eeg_values', 'outcome');
end

% Build the outcome array
%¨Params: 
% - recordings: times from the main input file
% - right: times only associated to the right hand event
% - left: times only associated to the left hand event
% r_or_l: identifying the hand targeted by the SOPHIA device (glove on the
%   right-'r' or left_'l' hand)
function outcome = buildOutcome(recordings, right, left, r_or_l)

    % Initialise the outcome array of length nbRecords in recordings with -1 values
    outcome = repmat(-1,1,length(recordings(:, 1)));

    switch r_or_l
        case 'r' %right hand targeted by the SOPHIA device
            outcome = buildHalf(recordings, right, 1, outcome);
            outcome = buildHalf(recordings, left, 0, outcome);
        case 'l' %left hand targeted by the SOPHIA device
            outcome = buildHalf(recordings, right, 0, outcome);
            outcome = buildHalf(recordings, left, 1, outcome);
        otherwise
            error('The third argument is either r (right hand) or l (left hand)\nYour parameter: %d', r_or_l);
    end
end

% Update the outcome by associating the right state (0 or 1) to the right
% cell, based on the recordings associated to one hand
% Params:
% - recordings: original times
% - hand: times associated to one hand
% - state: 1 (hand targeted by the SOPHIA device) or 0 (hand not targeted)
% - outcome: initial array of length nbRecords in recordings, filled with
% (-1) values
function outcome = buildHalf(recordings, hand, state, outcome)

    % Get the indices in recordings of every values shared with hand
    [~, indices] = ismember(hand, recordings);
    
    %Update the outcome array accordingly
    for i = indices
        outcome(i) = state;
    end
end