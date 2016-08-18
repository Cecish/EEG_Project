% Deciding which classifier to apply given an identifier
% Params:
%   - id: identifier of the classifier (1: kNN, 2: SVM, 3: MLP)
%   - final_mat_X: data used for the classification
%   - ex_events_Y: target field associated to each row of data
%   - k: number of nearest neighbours to consider for the kNN
% Return:  
%   - accuracy: accuracy score
%   - CVMdl: classifier model
function [accuracy, CVMdl] = classifiers( id, final_mat_X, ex_events_Y, k)
 
    switch id
        case 1 %kNN
            [accuracy, CVMdl] = kNN_func( final_mat_X, ex_events_Y, k );
        case 2 %SVM
            [accuracy, CVMdl] = SVM_func( final_mat_X, ex_events_Y );
%         case 3 %MLP
%             [accuracy, CVMdl] = MLP_func( final_mat_X, ex_events_Y );
        otherwise
            error('Wrong classifier identifier: %d', id);
    end
end


% Support Vector Machine
% Params: 
%   - final_mat_X: data used for the classification
%   - ex_events_Y: target field associated to each row of data
% Return: 
%   - accuracy: accuracy score
%   - CVSVMMdl: classifier model
function [accuracy, CVSVMMdl] = SVM_func( final_mat_X, ex_events_Y )

    %Cross Validation with the SVM classifier
    svmStruct = fitcsvm(final_mat_X, ex_events_Y);
    
    CVSVMMdl = crossval(svmStruct, 'KFold', 40);
    kloss = kfoldLoss(CVSVMMdl);
    accuracy = 100 - kloss*100;
    
    disp(['Accuracy = ', num2str(accuracy), '%']);
end


% K-Nearest Neighbours
% Params: 
%   - final_mat_X: data used for the classification
%   - ex_events_Y: target field associated to each row of data
%   - k: number of nearest neighbours to consider
% Return: 
%   - accuracy: accuracy score
%   - CVMdl: classifier model
function [accuracy, CVMdl] = kNN_func( final_mat_X, ex_events_Y, k )

    %Cross Validation with k-NN classifier
    Mdl = fitcknn(final_mat_X, ex_events_Y, 'NumNeighbors', k);
    CVMdl = crossval(Mdl, 'KFold', 40);
    kloss = kfoldLoss(CVMdl);
    accuracy = 100*(1 - kloss);
    
    disp(['Accuracy = ', num2str(accuracy), '%']);
end


% Multi-Layer Perceptron
% Params: 
%   - final_mat_X: data used for the classification
%   - ex_events_Y: target field associated to each row of data
% Return: accuracy score + trained ANN
function [accuracy, net] = MLP_func( final_mat_X, ex_events)
    
    k = 10; %k-folds for the cross-validation
    error = 0;
        
    for i = (1:k)
       display(['-------- Fold n°' , num2str(i), '--------'])
        
        %Define Multi-Layer Perceptron (creation + initialisation)
        testing_id = (i-1)*(40/k)+1:(i-1)*(40/k)+(40/k);
        training_id = setdiff(1: size(final_mat_X, 1),testing_id);
        net1 = newff(final_mat_X', ex_events, [10 10 10]);
        net1.divideFcn = 'divideind'; % Divide data by indices (i.e. not randomly)
        net1.divideParam.trainInd = training_id;
        net1.divideParam.valInd = [];
        net1.divideParam.testInd = testing_id;
        net1.trainParam.epochs = 200;

        % Train the Network
        [net, tr] = train(net1, final_mat_X', ex_events);

%       y = sim(net, final_mat_X');
%       min(tr.tperf)
        error = error + min(tr.tperf);
    end
    accuracy = (1 - error/40)*100;
    
    % Test the Network
%       y = sim(net, final_mat_X');
%      accuracy = 100 - immse(ex_events, y);
    disp(['Accuracy = ', num2str(accuracy), '%']);
end


% Applying the kNN on a traning set and getting an accuracy score
% Params:
%   - final_mat_X: eeg recordings, 2D matrix (nb_records x channels)
%   - ex_events_Y: hand state associated to each recordings
%   - nb_trials : number of trials for building the training set from the
%   whole data matrix
%   - tot_trials: number of trials when acquiring the data
%   - k: number of neighbours to consider
% Retun : accuracy score and predictions made
% function [predictions, accuracy] = kNN( final_mat_X, ex_events_Y, nb_trials, tot_trials, k )
% 
%     % #### 1: Min-Max normalisation
%     width = size(final_mat_X, 2); %14, 24
%     height = size(final_mat_X, 1);   %3280, 560
%     manipFuns = dataManipFunctions;
%     minmax_X = manipFuns.MinMaxNorm(final_mat_X, height, width);
%  
%     % #### 2: Randomly shuffle trials in final_mat_X and ex_events_Y
%     %[X, Y] = randomShuffle(minmax_X, ex_events_Y, tot_trials);
%     
%     % #### 3: Split data into test train/train dataset
%     [training_dataset, testing_dataset, training_Y, testing_Y] = ...
%         manipFuns.splitXY(minmax_X, ex_events_Y, nb_trials, tot_trials);
%     %[training_dataset, testing_dataset, training_Y, testing_Y] = ...
%     %    manipFuns.splitXY(X, Y, nb_trials, tot_trials);
% 
%     dlmwrite('enormetest2.txt', [final_mat_X ex_events_Y'], 'delimiter', ',');
%     
%     % #### 4: Apply kNN on test set with train set as reference
%     predictions = aux_knn(training_dataset, testing_dataset, k, training_Y);
%     accuracy = manipFuns.calculateAccuracy(predictions, testing_Y, ...
%         tot_trials - nb_trials);
%     disp(['Accuracy with ', num2str(k), '-NN = ', num2str(accuracy), '%']);
%     
%     % Proportion
%     %zero_y = sum(ex_events_Y ~= 0)/height*100;
%     zero_y = sum(testing_Y ~= 0)/(tot_trials - nb_trials)*100;
%     %disp(['1 = ', num2str(zero_y), '% and 0 = ', num2str(100-zero_y), '%']);
% end


% Randomly shuffle the dataset by trials
% Params:
%   - final_mat_X: matrix to shuffle
%   - ex_events: associated events to shuffle the same way
%   - tot_trials: total number of trials
% Return:
%   X: data matrix shuffled
%   Y: events array shuffled accordingle
% function [X, Y] = randomShuffle(final_mat_X, ex_events, tot_trials)
%     
%     %Randomly decide on the new trials order
%     ordered_trials = [1:tot_trials];
%     shuffled_trials = ordered_trials(randperm(length(ordered_trials)));
%     
%     %Update (final_mat_X, ex_events) accordingly in (X, Y)
%     [X, Y] = updateShuffle(final_mat_X, ex_events, tot_trials, shuffled_trials);
% end


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
% function [X, Y] = updateShuffle(mat_X, events, tot_trials, shuffled_trials)
% 
%     my_length = 0;
%     nb_records_per_trial = size(mat_X, 1)/tot_trials;
%     X = [];
%     Y = [];
%     
%     while ~isequal(my_length, tot_trials)
%         pos_i = ((shuffled_trials(my_length+1)-1) * nb_records_per_trial) + 1;
%         
%         temp = mat_X((pos_i : pos_i + nb_records_per_trial - 1), :);
%         X = vertcat(X, temp);
%         
%         tempY = events(pos_i : pos_i + nb_records_per_trial - 1);
%         Y = vertcat(Y, tempY');
%         
%         my_length = my_length + 1;
%     end
% end


% Extract the k nearest neighbours to a test instance
% Params: 
%   - training_set
%   - test_instance: test data point
%   - k: number of neighbours to consider 
% Return: k rows() in the training set that are(is) the closest to the test
% instance according to the euclidian distance
% function neighbours = getNeighbours(training_set, test_instance, k)
% 
%     for x = (1:size(training_set, 1))
% 		distances(x, :) = [sqrt(sum((test_instance - training_set(x, :))...
%             .^ 2)) training_set(x, :)];
%     end
%     
% 	distances = sortrows(distances,1);
% 	neighbours = [];
% 	for x = (1:k)
% 		neighbours(x, :) = distances(x, 2:size(distances, 2));
%     end
% end


% Calculate the accuracy score
% Params: 
%   - training_dataset
%   - testing_dataset
%   - k: number of nearest neighbours considered
%   - training_Y: events associated to the training data
% Return: class predictions
% function predictions = aux_knn(training_dataset, testing_dataset, k, training_Y)
%     temp = 0.0;
%         
%     for i = (1: size(testing_dataset, 1))
%         neighbours = getNeighbours(training_dataset, testing_dataset(i, :), k);
%         
%         % Get the index of each neighbours extracted in the training dataset
%         for j = (1:k)
%             idx(j) = find(ismember(training_dataset, neighbours(j, :)),1);
%         end
%         
%         % Get the number of 1- events associated to the neighbours
%         nb_1 = sum(training_Y(idx) ~= 0);
%         % Decide of the state of the current test instance by taking the
%         % mean of each neighbours' events
%         if ((nb_1/length(idx)*100) > 50)
%             state = 1.0;
%         else
%             state = 0.0;
%         end
%         
%         predictions(i) = state;
%     end
% end


% function [accuracy, net] = CV_test( dataset, net, k)
%     count = 0;
%     
%     %1. Divide the data into k non-overlapping folds
%     %class 0
%     struct_folds1 = selectClass(dataset, k/2, 0);
%     %lass 1
%     struct_folds2 = selectClass(dataset, k/2, 1);
%     %Random order for the k folds
%     ordering = randperm(size(dataset, 1)); %4
%     temp = [struct_folds1'; struct_folds2'];
%     %struct_folds = shuffleFolds(temp, ordering, k);
%     % Transform temp struct in matrix because so much easy to manipulate
%     idx = 1;
%     for p = (1: k) %10 folds
%         for r = (1: 4)
%            temp2(idx, :) = temp{p}(r, :);
%            idx = idx + 1;
%         end
%     end
% 
%     temp_mat = temp2(ordering, :);
%     %Deconstruct in struct folds
%     for pp = (1: k)
%         struct_folds{pp} = temp_mat((pp-1)*4+1: (pp-1)*4 + 4, :);
%     end
%     
%     % Create a Pattern Recognition Network
%     hiddenLayerSize = 10;
%     net = patternnet(hiddenLayerSize);
%     net.numLayers = 3;
%     net.trainParam.showWindow = false;
%     net.trainParam.showCommandLine = false; 
%     net
%     
%     %For each fold
%     for i = (1: k)
%         training_set = buildTrainingSet(struct_folds, i);
%         idx = 1;
%         for p = (1: 9) %10 folds
%             for r = (1: 4)
%                training_set2(idx, :) = training_set{p}(r, :);
%                idx = idx + 1;
%             end
%         end
% 
%         % Setup Division of Data for Training, Testing
%         net.divideParam.trainRatio = 70/100;
%         %net.divideParam.valRatio = 15/100;
%         net.divideParam.testRatio = 30/100;
% 
%         % Train the Network
%         [net,tr] = train(net, training_set2(:, 1:size(training_set2, 2)-1)', training_set2(:, size(training_set2, 2))');
%         
%         % Test the classier on all the examples in Fold i
%         y = net(struct_folds{i}(:, 1: size(struct_folds{i}, 2)-1)');
%         %predictions = decisionOutcome(y, length(y));
%               
%         % Number of records wrongly classified in Fold i
%         mse = 0;
%         outcome = struct_folds{i}(:, size(struct_folds{i}, 2));
%         for z = 1: length(y)
%             mse = mse + (y(z) - struct_folds{i}(z, size(struct_folds{i}, 2)))^2;
%         end
%         res = mse/length(y);
%     end
%     accuracy = res/40;
% end

% function training_set = buildTrainingSet(struct_folds, id)
%     training_set = [];
%     
%     for i = (1: 10)
%         if ~isequal(i, id)
%            training_set = [training_set; struct_folds(i)]; 
%         end
%     end
% end

% function struct_folds = shuffleFolds(my_struct, ordering, k)
%     count = 0;
%     idx = 1;
%     
%     for i = (1: length(ordering))
%        count = count +1;
%        if isequal(count,  5)
%            idx = idx + 1;
%            count = 1;
%        end
%        disp([num2str(idx), ' - ', num2str(count), ' - ', num2str(round(ordering(i)/(40/k))+1), ' - ', num2str(mod(ordering(i), 40/k)+1)]);
%        struct_folds{idx}(count, :) = my_struct{round(ordering(i)/(40/k))+1}(mod(ordering(i), 40/k)+1, :); 
%     end
% end


% Select a specific number of rows of a specific class
% Params:
%   - dataset: (data + target field as last column)
%   - k: number of folds/2
%   - class_id: either 1 or 0 (correct hand or not)
% Return:
%   - struct_folds: k-folds structure
% function struct_folds = selectClass(dataset, k, class_id)
% 
%     % Select rows with the right target field
%     indices = dataset(:, size(dataset, 2)) == class_id;
%     sub_dataset = dataset(indices, :);
% 
%     % Random orders
%     orders = randperm(size(sub_dataset, 1));
% 
%     % Build sub matrix of selected rows & Update initial dataset
%     for i = (1: k)
%         my_temp = [];
%         for z = ((i-1)*(4)+1: (i-1)*(4)+4)
%             my_temp(z - (i-1)*(4), :) = sub_dataset(orders(z), :);
%         end
%         struct_folds{i} = my_temp;
%     end
% end
