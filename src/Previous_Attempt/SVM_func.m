function [predictions, accuracy] = SVM_func( final_mat_X, ex_events_Y, ...
    nb_trials, tot_trials )

    % #### 1: Min-Max normalisation
    width = size(final_mat_X, 2); %14, 24
    height = size(final_mat_X, 1);   %3280, 560
    manipFuns = dataManipFunctions; 
    
    minmax_X = manipFuns.MinMaxNorm(final_mat_X, height, width);
    
    % #### 3: Split data into test train/train dataset
    [training_dataset, testing_dataset, training_Y, testing_Y] = ...
        manipFuns.splitXY(minmax_X, ex_events_Y, nb_trials, tot_trials);
    
    %Train SVM
    svmStruct = svmtrain(training_dataset, training_Y, 'ShowPlot',false);
    
    %Test SVM
    predictions = svmclassify(svmStruct, testing_dataset, 'ShowPlot', false);
    
    accuracy = manipFuns.calculateAccuracy(predictions, testing_Y, ...
        tot_trials-nb_trials);
    disp(['Accuracy = ', num2str(accuracy), '%']);
end


