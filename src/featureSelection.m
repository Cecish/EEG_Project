function res = featureSelection( mat, pop_size, nb_features, ...
        nb_iterations, k, nb_trials, nb_trials_training, events,...
        crossover_rate, mutation_rate)
    best_features = [];
    it = 0;
    best_fitness = -Inf;
    fitness_evolution = [best_fitness];
        
    % Randomness using seed
    %rand('state',sum(100*clock));
    
    % Population initialisation
    pop = popInit(pop_size, nb_features);
        
    %Loop until the termination condition is reached
    while ((it < nb_iterations) & (best_fitness < 99.99))
        disp(['#### EPOCH n° ', num2str(it)]);

        % ####1: Evaluate the fitness of each chromosome in the population
        fitness_array = evalPop(mat, pop, pop_size, nb_features, events, ...
            nb_trials_training, nb_trials, k);

        % Sort the population by fitness and record the best fitness of the generation
        [pop, fitness_evolution, best_fitness, best_features] = sortRecord(pop, fitness_array, ...
            fitness_evolution, best_fitness, best_features);
        %size(best_features)
        % ####2: Create a new population
        pop = reproduction(pop, pop_size, crossover_rate, mutation_rate, ...
            fitness_array, nb_features);
        
        it = it + 1; % Increment
    end
    
    %Save the fitness evolution in a file
    dlmwrite('fitness_evol.csv', fitness_evolution', 'delimiter', ';');
    
    disp(['Best features: ', mat2str(best_features)]);
    
    %Build matrix with the best features
    res = buildSubMat(mat, best_features, nb_features);
end


% Randomly create a chromosome
% Param: nb_features: number of genes for the chromosome
% Return: new chromosome
function pop = createIndividual(nb_features)
    for i = (1: nb_features)
        if (rand()<0.5)
            pop(i) = 0;
        else
            pop(i) = 1;
        end
    end
end


% Given a features binary array, a sub matrix of mat is built
% Params: 
%   - mat: original matrix with all the features
%   - individual: features binary array
%   - length_id: length of the features array
% Return: sub matrix of mat
function new_mat = buildSubMat(mat, individual, length_id)
    count = 1;

    for i = (1 : length_id)
        if isequal(individual(i), 1) %Features taken for new_mat
            new_mat(:, count) = mat(:, i);
            count = count + 1;
        end
    end
end


% Selection with the tournament selection
% Params: 
%   - pop: population of chromosmes
%   - pop_size: size of the population
%   - nb: number of chrosomes to pre-select (4 is a good number for a 
%        population of size 200)
%   - fitness_array
% Return: the two parents (chromosomes) selected
function [parent1, parent2] = selection(pop, pop_size, nb, fitness_array)
    
    %Select nb random individuals in the population
    random2 = Inf;
    for l = (1:nb)
        random = round(rand()*(pop_size-1) + 1);
        while isequal(random, random2)
            random = round(rand()*(pop_size-1) + 1);
        end
        random2 = random;

        temp(l, :) = [pop(random, :) fitness_array(random)];
    end
    
    %Sort selected individuals by fitness
    [~, idx] = sort(temp(:, size(temp, 2)), 'descend');
    temp = temp(idx, :);
    
    % Select the best two chromosomes out of the pre-selected subset
    parent1 = temp(1, 1:length(temp)-1);
    parent2 = temp(2, 1:length(temp)-1);
end


% One point crossover operator
% Params: 
%   - parent1
%   - parent2
%   - length: length of the chromosomes
% Return: 
%   - child1: of parent1 and parent2 (inheritance of their genes)
%   - child2: of parent1 and parent2 (inheritance of their genes)
function [child1, child2] = onePointCrossover(parent1, parent2, length)

    % Defining the crossover point randomly
    rand_pos = round((length) *rand());
    
    for i = (1:length)
        %Before crossover point, child1 inherits from parent1's genes and
        %child2 inherits from parent2's genes
       if (i < rand_pos) 
           child1(i) = parent1(i);
           child2(i) = parent2(i);
        
       %After crossover point, child1 inherits from parent2's genes and 
       %child2 inherits from parent2's genes
       else
           child1(i) = parent2(i);
           child2(i) = parent1(i);
       end
    end
end


% Two points crossover operator
% Params: 
%   - parent1
%   - parent2
%   - length: length of the chromosomes
% Return: 
%   - child1: of parent1 and parent2 (inheritance of their genes)
%   - child2: of parent1 and parent2 (inheritance of their genes)
function [child1, child2] = twoPointCrossover(parent1, parent2, length)
    %Randomly defining the two points of crossover
    rand_pos1 = round((length-1) *rand()) + 1;
    rand_pos2 = rand_pos1;
    
    %Ensuring that the two points are distinct
    while isequal(rand_pos1, rand_pos2)
        rand_pos2 = round((length-1) *rand()) + 1;
    end
    
    % Genes inheritance process
    for i = (1:length)
       if ((i < min(rand_pos1, rand_pos2)) || (i> max(rand_pos1, rand_pos2))) 
           child1(i) = parent1(i);
           child2(i) = parent2(i);
       else
           child1(i) = parent2(i);
           child2(i) = parent1(i);
       end
    end
end


% Mutation operator. The chromosome is a binary vector, so the mutation
% consists in flipping specific genes (1->0 or 0->1)
% Params: 
%   - chromosome: to be mutated
%   - rate: for determining how much genes to mutate
%   - length: length of the chromosome
% Return: the mutated chromosome
function chromosome = mutation(chromosome, rate, length)
    nb_to_mutate = (length * rate); %Number of genes to mutate
    
    for i = (1: nb_to_mutate)
        %Select genes to mutate randomly
        temp = rand() * (length-1) + 1; 
        rand_pos_m = round(temp);
        
        % Flips gene value from 1 to 0
        if isequal(chromosome(rand_pos_m),1)
            chromosome(rand_pos_m) = 0;
        % Flips gene value from 0 to 1
        else    
            chromosome(rand_pos_m) = 1;
        end
    end
end


% Initial population generation
% Params: 
%   - pop_size: size of the population to create
%   - nb_features: number of genes for each individual
% Return : the new population
function pop = popInit(pop_size, nb_features)
    for i = (1: pop_size)
        pop(i, :) = createIndividual(nb_features);
    end
end


% Evaluate the current population
% Params: 
%   - mat: matrix of features
%   - pop: population of chromosomes
%   - pop_size: size of the population
%   - nb_features: number of features for each chromosomes
%   - events: events associated to each row of mat
%   - nb_trials_training: number of trials that constitute the "training" 
%       dataset for the k-NN
%   - nb_trials: total number of trials
%   - k: number of neighbours to consider for the k-NN
% Return: a fitness array where each value corresponds to the fitness of
% each chromosome of the population
function fitness_array = evalPop(mat, pop, pop_size, nb_features, events, ...
    nb_trials_training, nb_trials, k)
    
    for i = (1: pop_size)
        %Build new mat according to features selected
        new_mat = buildSubMat(mat, pop(i, :), nb_features);
            
        %[~, accuracy, ~] = kNN(new_mat, events, nb_trials_training, nb_trials, k);
        accuracy = SVM_func(new_mat, events, nb_trials_training, nb_trials);
        fitness_array(i) = accuracy;
    end
end


% Sort the population by fitness (descending order) and keep a record of
% the best fitness achieved for this generation
% Params: 
%   - pop: population of chromosomes
%   - fitness_evolution: array of best fitness per generation
%   - best fitness of the generation
%   - best_features: best chromosome of the population
% Return: 
%   - sorted population
%   - (updated) best fitness
%   - updated array of best fitness per generation
%   - best_features: best individual of the population
function [pop, fitness_evolution, best_fitness, best_features] = sortRecord(pop, ...
    fitness_array, fitness_evolution, best_fitness, best_features)
    
    pop = [pop fitness_array'];
    
    %sortrows(pop, size(pop, 2), 'descend');
    [~, idx] = sort(pop(:, size(pop, 2)), 'descend');
    pop = pop(idx, :);
        
    if (pop(1, size(pop, 2)) > best_fitness)
        best_fitness = pop(1, size(pop, 2));
        best_features =  pop(1, 1:size(pop, 2)-1);
    end
    
    fitness_evolution = [fitness_evolution best_fitness];
    pop = pop(:, 1 : size(pop, 2)-1);
end
        

% Population reproduction
% Params: 
%   - pop: population of chromosomes
%   - pop_size: size of the population
%   - crossover_rate
%   - mutation_rate 
%   - fitness_array
%   - nb_features
% Return: new population
function pop = reproduction(pop, pop_size, crossover_rate, mutation_rate, ...
            fitness_array, nb_features)
    
    % N best individuals to directly copy/paste in the new generation    
    delimiter_unchanged = round(pop_size * (1-crossover_rate));
    
    % For each individual
    for i = (1: pop_size)
        %Directly copy/paste
        if (i<= delimiter_unchanged)
            new_pop(i, :) = pop(i, :);

        else
            % 2.1 Selection: tournament selection
            [parent1, parent2] = selection(pop, pop_size, 4, fitness_array);
        
            % 2.2 Crossover (one point)
            %[child1, child2] = onePointCrossover(parent1, parent2, nb_features);
            [child1, child2] = twoPointCrossover(parent1, parent2, nb_features);

            % 2.3 Mutation
            child1 = mutation(child1, mutation_rate, nb_features);
            child2 = mutation(child2, mutation_rate, nb_features);
               
            new_pop(i, :) = child1;
            new_pop(i+1, :) = child2;
        end
    end
    
    % Ensuring that the sizes of the new_population and the old population
    % are the same
    while (size(new_pop, 1) > pop_size)
        new_pop = new_pop(1:size(new_pop, 1)-1, :); 
    end
    
    %Update the population
    pop = new_pop;    
end
   
