nb_runs = 2;
% 
% save('input_params', readParameterFile());

% #### 0. Preprocessing ####
% 0.1 Filtering
accuracies = zeros(1, nb_runs);
for i = (1:nb_runs)
%     clear all;
%     load input_params;
    
    display([' #### Iteration n° ', num2str(i), ' ####']);
    accuracies(i) = main();
end

accuracies
mean(accuracies)
var(accuracies)
std(accuracies)

dlmwrite('accuracies.csv', accuracies, 'delimiter', ';');