% Functions are fields of a struct (holder)
% Return function handles to local functions
function funs = dataManipFunctions
    funs.MinMaxNorm = @MinMaxNorm; %Min-max normalisation
    funs.splitXY = @splitXY; %split data matrix into train and test sub-matrices
    funs.calculateAccuracy = @calculateAccuracy; %Accuracy score calculation
    funs.featureSelection = @featureSelection; %Choice of the method to use for feature selection
    funs.analysisSelectedFeatures = @analysisSelectedFeatures; %Analysis of the features selected
end


% Min max normalisation of a 2D matrix
% Params: 
%   - dataset: 2D matrix
%   - height: of the matrix (number of rows)
%   - width: of the matrix (number of columns)
% Return: the coresponding minmax normalised matrix
function res_matrix = MinMaxNorm(dataset, height, width)

    res_matrix = zeros(height, width); %Initialisation

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


% Split the data matrix and events into train and test subsets
% Ratio: nb_trials/(total_trials _ nb_trials) for train/test
% Params: 
%   - X: data matrix
%   - Y: corresponding array of events
%   - nb_trials: number of trials in the train subset
%   - tot_trials: total number of trials
% Return: 
%   - training_dataset
%   - testing_dataset
%   - training_Y: extracted from the events array
%   - testing_Y: extracted from the events array
function [training_dataset, testing_dataset, training_Y, testing_Y] = ...
    splitXY(X, Y, nb_trials, tot_trials)

    nb_rows_per_trial = size(X, 1) / tot_trials;
    
    training_dataset(1:nb_trials * nb_rows_per_trial, :) = X(1:nb_trials *...
        nb_rows_per_trial, :);
    testing_dataset(1:(size(X, 1)-(nb_trials * nb_rows_per_trial)), :) = ...
        X((nb_trials * nb_rows_per_trial)+1:size(X, 1), :);
    training_Y(1:nb_trials * nb_rows_per_trial) = Y(1:nb_trials * ...
        nb_rows_per_trial);
    testing_Y(1:(size(X, 1)-(nb_trials * nb_rows_per_trial))) = Y((nb_trials...
        * nb_rows_per_trial)+1:size(X, 1));
end


% Calculation of the accuracy score
% Params: 
%   - predictions: classes predicted by the classifier
%   - observed: correct classes
%   - length: number of classes predicted and observed
% Return: the accuracy score given the predictions and the expected
% observations
function accuracy = calculateAccuracy(predictions, observed, length)
    temp = 0.0;
    
    for i = (1: length)
       if isequal(predictions(i), observed(i))
          temp = temp + 1;
       end
    end
    
    accuracy = (temp / length) * 100;
end


% Analysis of the best features selected per channel by plotting a stack histogram
% Params: 
%   - features: best features selected with a genetic algorithm
%   - nb_channels: number of channels considered for the analysis
%   - ratio_features: number of features extracted in each feature category
%   (DWT, AR, PSD, Time domain features)
function analysisSelectedFeatures(features, nb_channels, ratio_features)
    % Number of different types of features extracted
    nb_cat_features = 2; %wavelet and AR features
    nb_features_per_channel = length(features)/nb_channels; % 100%
  
    %Initialisation of the matrix used for plotting the histogram
    mat_features = zeros(nb_channels, nb_cat_features);
    
    % For each channel
    for i = (1: nb_channels)
        %Number of wavelet features selected
        mat_features(i, 1) = sum(features((i-1)*nb_features_per_channel+1:(i-1)...
            *nb_features_per_channel+ratio_features(1)) == 1);
        %Number of AR features selected
        mat_features(i, 2) = sum(features((i-1)*nb_features_per_channel+...
            ratio_features(1)+1:(i-1)*nb_features_per_channel+ratio_features(1)...
            +ratio_features(2)) == 1);
        %Number of PSD features selected
        mat_features(i, 3) = sum(features((i-1)*nb_features_per_channel+...
            ratio_features(1)+ratio_features(2)+1:(i-1)*nb_features_per_channel...
            +ratio_features(1)+ratio_features(2)+ratio_features(3)) == 1);
        %Number of time domain (TD) features selected
        mat_features(i, 4) = sum(features((i-1)*nb_features_per_channel+...
            ratio_features(1)+ratio_features(2)+ratio_features(3)+1:(i-1)*...
            nb_features_per_channel+ratio_features(1)+ratio_features(2)+...
            ratio_features(3)+ratio_features(4)) == 1);
    end
    
    % The 4 types of selected features displayed in the stacked bar chart
    % are represented by a percentage. Consequently the maximum height that
    % can possibly be reached in the stacked bar chart is 400 = 4*1000%
    mat_features_bis(:, 1) = mat_features(:, 1)*100/ratio_features(1);
    mat_features_bis(:, 2) = mat_features(:, 2)*100/ratio_features(2);
    mat_features_bis(:, 3) = mat_features(:, 3)*100/ratio_features(3);
    mat_features_bis(:, 4) = mat_features(:, 4)*100/ratio_features(4);
    
    % Plotting the stacked bar chart
    figure('Name', 'Features distribution')
    bar_plot = bar(mat_features_bis, 'stacked');
    title('Features distribution')
    xlabel('Channels')
    ylabel('Cummulative percentage')
    grid on
    legend(bar_plot, {'Discrete Wavelet Transform (DWT) features', ...
        'Auto-regressive (AR) model features', ...
        'Power Spectrum Analysis (PSD) features', 'Time Domain (TD) features'},...
        'Location','Best','FontSize',8);
end