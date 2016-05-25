% Applying the kNN on a traning set and getting an accuracy score
% Params:
%   - final_mat_X: eeg recordings, 2D matrix (nb_records x channels)
%   - ex_events_Y: hand state associated to each recordings
%   - nb_trials : number of trials for building the training set from the
%   whole data matrix
%   - tot_trials: number of trials when acquiring the data
%   - k: number of neighbours to consider
% Retun : accuracy score
function [predictions, accuracy, observed] = kNN( final_mat_X, ex_events_Y, nb_trials, tot_trials, k )

    % #### 1: Min-Max normalisation
    width = size(final_mat_X, 2); %14, 24
    height = size(final_mat_X, 1);   %3280, 560
    minmax_X = MinMaxNorm(final_mat_X, height, width);
 
    % #### 2: Randomly shuffle trials in final_mat_X and ex_events_Y
    %[X, Y] = randomShuffle(minmax_X, ex_events_Y, tot_trials);
    
    % #### 3: Split data into test train/train dataset
    [training_dataset, testing_dataset, training_Y, testing_Y] = splitXY(...
        minmax_X, ex_events_Y, nb_trials, tot_trials);
    %[training_dataset, testing_dataset, training_Y, testing_Y] = splitXY(...
    %    X, Y, nb_trials, tot_trials);

    dlmwrite('enormetest2.txt', [final_mat_X ex_events_Y'], 'delimiter', ',');
    
    % #### 4: Apply kNN on test set with train set as reference
    [predictions, accuracy, observed] = aux_knn(training_dataset, testing_dataset, k, training_Y, testing_Y);
    disp(['Accuracy with ', num2str(k), '-NN = ', num2str(accuracy), '%']);
    
    % Proportion
    %zero_y = sum(ex_events_Y ~= 0)/height*100;
    zero_y = sum(testing_Y ~= 0)/(tot_trials - nb_trials)*100;
    %disp(['1 = ', num2str(zero_y), '% and 0 = ', num2str(100-zero_y), '%']);
end

% Min max normalisation of a 2D matrix
% Params: 
%   - dataset: 2D matrix
%   - height: of the matrix
%   - width: of the matrix
% Return: the coresponding minmax normalised matrix
function res_matrix = MinMaxNorm(dataset, height, width)

    %Finding the min and the max valus of each columns
    minmax_arrays = findMinMax(dataset, height, width);

    for j = (1: width)
        for i = (1 : height)
            res_matrix(i, j) = (dataset(i, j)-minmax_arrays(1, j))/...
                (minmax_arrays(2, j)-minmax_arrays(1, j));
        end
    end
end

% Find min max arrays for each column of the data matrix
% Params: 
%   - dataset: matrix
%   - height: of the matrix
%   - width: of the matrix
% Return: A 2 x nb_col matrix where the first row stores the min of each
% column and the second row stores the max of each column
function minMax_arrays = findMinMax(dataset, height, width)
    
    %Min
    for j = (1: width)
        min = Inf;
        for i = (1:height)
            if (dataset(i, j) < min)
                min = dataset(i, j);
            end
        end
        
        %disp(['Minimum of column', num2str(j), ' = ', num2str(min)]);
        minMax_arrays(1, j) = min;
    end

    %Max
    for i = (1 : width)
        max = -Inf;
        for j = (1: height)
            if ( dataset(j, i) > max)
                max = dataset(j, i);
            end
        end
        
        %disp(['Maximum of column', num2str(i), ' = ', num2str(max)]);
        minMax_arrays(2, i) = max;
    end
end

% Randomly shuffle the dataset by trials
% Params:
%   - final_mat_X: matrix to shuffle
%   - ex_events: associated events to shuffle the same way
%   - tot_trials: total number of trials
% Return:
%   X: data matrix shuffled
%   Y: events array shuffled accordingle
function [X, Y] = randomShuffle(final_mat_X, ex_events, tot_trials)
    
    %Randomly decide on the new trials order
    ordered_trials = [1:tot_trials];
    shuffled_trials = ordered_trials(randperm(length(ordered_trials)));
    
    %Update (final_mat_X, ex_events) accordingly in (X, Y)
    [X, Y] = updateShuffle(final_mat_X, ex_events, tot_trials, shuffled_trials);
end


% Auxiliary function to randomShuffle that given a new order od trials will
% shuffle both the data matrix and the events array
% Params: 
%   - mat_X: data matrix
%   - events: events array
%   - tot_trials: total number of trials
%   - shuffled_trials: new order of trials generated randomly
% Return: 
%   - X: shuffled data matrix
%   - Y: suffled events array
function [X, Y] = updateShuffle(mat_X, events, tot_trials, shuffled_trials)

    my_length = 0;
    nb_records_per_trial = size(mat_X, 1)/tot_trials;
    X = [];
    Y = [];
    
    while ~isequal(my_length, tot_trials)
        pos_i = ((shuffled_trials(my_length+1)-1) * nb_records_per_trial) + 1;
        
        temp = mat_X((pos_i : pos_i + nb_records_per_trial - 1), :);
        X = vertcat(X, temp);
        
        tempY = events(pos_i : pos_i + nb_records_per_trial - 1);
        Y = vertcat(Y, tempY');
        
        my_length = my_length + 1;
    end
end


% Split the data matrix and events into train and test subsets
% Ratio: nb_trials/(total_trials _ nb_trials) for train/test
% Params: 
%   - X: data matrix
%   - Y: corresponding events array
%   - nb_trials: number of trials in the train subset
%   - tot_trials: total number of trials
% Return: 
%   - training_dataset
%   - testing_dataset
%   - training_Y: extracted from the events array
%   - testing_Y: extracted from the events array
function [training_dataset, testing_dataset, training_Y, testing_Y] = splitXY(X, Y, nb_trials, tot_trials)

    nb_rows_per_trial = size(X, 1) / tot_trials;
    
    training_dataset(1:nb_trials * nb_rows_per_trial, :) = X(1:nb_trials * nb_rows_per_trial, :);
    testing_dataset(1:(size(X, 1)-(nb_trials * nb_rows_per_trial)), :) = X((nb_trials * nb_rows_per_trial)+1:size(X, 1), :);
    training_Y(1:nb_trials * nb_rows_per_trial) = Y(1:nb_trials * nb_rows_per_trial);
    testing_Y(1:(size(X, 1)-(nb_trials * nb_rows_per_trial))) = Y((nb_trials * nb_rows_per_trial)+1:size(X, 1));
end


% Extract the k nearest neighbours to a test instance
% Params: 
%   - training_set
%   - test_instance: test data point
%   - k: number of neighbours to consider 
% Return: k rows() in the training set that are(is) the closest to the test
% instance according to the euclidian distance
function neighbours = getNeighbours(training_set, test_instance, k)

    for x = (1:size(training_set, 1))
		distances(x, :) = [sqrt(sum((test_instance - training_set(x, :))...
            .^ 2)) training_set(x, :)];
    end
    
	distances = sortrows(distances,1);
	neighbours = [];
	for x = (1:k)
		neighbours(x, :) = distances(x, 2:size(distances, 2));
    end
end

% Calculate the accuracy score
% Params: 
%   - training_dataset
%   - testing_dataset
%   - k: number of nearest neighbours considered
%   - training_Y: events associated to the training data
%   - testing_Y: events associated to the testing data
% Return: the accuracy score
function [predictions, accuracy, testing_Y] = aux_knn(training_dataset, testing_dataset, k, training_Y, testing_Y)
    temp = 0.0;
        
    for i = (1: size(testing_dataset, 1))
        neighbours = getNeighbours(training_dataset, testing_dataset(i, :), k);
        
        % Get the index of each neighbours extracted in the training dataset
        for j = (1:k)
            idx(j) = find(ismember(training_dataset, neighbours(j, :)),1);
        end
        
        % Get the number of 1- events associated to the neighbours
        nb_1 = sum(training_Y(idx) ~= 0);
        % Decide of the state of the current test instance by taking the
        % mean of each neighbours' events
        if ((nb_1/length(idx)*100) > 50)
            state = 1.0;
        else
            state = 0.0;
        end
        
        predictions(i) = state;
        
        %If guess is correct, update the accuracy
        if isequal(state, testing_Y(i))
            temp = temp + 1;
        end
    end
    
    %accuracy = temp / length(training_Y)*100;
    accuracy = temp / length(testing_Y)*100;
end
