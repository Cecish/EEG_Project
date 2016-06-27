% Feature selection
% Params: 
%   - mat: initial matrix with all the features previously extracted
%   - nb_features: total number of features
%   - k: number of nearest neighbours to consider with the k-NN classifier
%   - events: known outcome associated to each trial
%   - classifier: id of the chosen classifier 
% Return: 
%   - res: sub matrix (with a reduced number of features) that led to the
%     best accuracy score
%   - best_features: subset of features that led to the best accuracy score
%   - best_net: trained MLP that led to the best accuracy score (when the MLP 
%     classifier is chosen)
function [res, best_features] = featureSelection( mat, nb_features, k,...
    events, classifier)
    
    % GA parameters
    pop_size = 20;
    max_generations = 5;
    nb_elites = 2;
    crossover_rate = 0.8;
    mutation_rate = 0.1;
    
    % Parameters initialisation
    best_features = [];
    it = 0;
    best_fitness = -Inf;
    fitness_evolution = best_fitness;
        
    % Randomness using seed
    %rand('state',sum(100*clock));
    
    % Population initialisation
    pop = popInit(pop_size, nb_features);
        
    %Loop until the termination condition is reached
    while ((it < max_generations) && (best_fitness < 100.00))
        disp(['#### EPOCH n° ', num2str(it)]);

        % ####1: Evaluate the fitness of each chromosome in the population
        fitness_array = evalPop(mat, pop, pop_size, ...
            nb_features, events, k, classifier);

        % Sort the population by fitness and record the best fitness of the generation
        [pop, fitness_evolution, best_fitness, best_features, fitness_array] = ...
            sortRecord(pop, fitness_array, fitness_evolution, ...
            best_fitness, best_features);

        % ####2: Create a new population
        pop = reproduction(pop, pop_size, crossover_rate, mutation_rate, ...
            fitness_array, nb_features, nb_elites);
        
        it = it + 1; % Increment
        display(['Iteration: ',num2str(it),': Best fitness = ', ...
            num2str(fitness_evolution(it))]);
    end
    
    %Save the fitness evolution in a file
    dlmwrite('fitness_evol.csv', fitness_evolution', 'delimiter', ';');
    
    disp(['Best features: ', mat2str(best_features)]);
    
    %Build matrix with the best features
    res = buildSubMat(mat, best_features, nb_features);
    
    figure('Name', 'Fitness (accuracy) evolution over the generations')
    plot(1:it, fitness_evolution(2:end))
    title('Best fitness evolution over the generations')
    xlabel('Generations')
    ylabel('Fitness, Accuracy')
    
%     figure('Name', 'Features distribution')
%     bar(best_features, 1.2) %Some space added to see the top of the bars (aesthetic purpose)
%     title('Features distribution')
%     xlabel('Features')
%     ylabel('Apparition')
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
function [parent1, parent2] = tournamentSelection(pop, pop_size, nb, fitness_array)
    
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


% Uniform crossover. A binary crossove rmask is created. The first child
% inherits the genes from the first parent every time the mask components are
% true, and inherits genes from the second parent where the mask's 
% componenents are false
% Params: 
%   - parent1
%   - parent2
% Return: 
%   - child1: of parent1 and parent2 (inheritance of their genes)
%   - child2: of parent1 and parent2 (inheritance of their genes)
function [child1, child2] = uniformCrossover(parent1, parent2, length)
    % Creation of the mask
    mask = randi([0 1], 1, length);
    child1=mask.*parent1+(1-mask).*parent2;
    child2=mask.*parent2+(1-mask).*parent1;
end


% Mutation operator. The chromosome is a binary vector, so the mutation
% consists in flipping specific genes (1->0 or 0->1)
% Params: 
%   - chromosome: to be mutated
%   - rate: for determining how much genes to mutate
%   - length: length of the chromosome
% Return: the mutated chromosome
function chromosome = mutation(chromosome, rate, length)
    
    for i = (1: length)
        if(rand() <= rate)
           chromosome(i) = 1-chromosome(i); 
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
%   - k: number of neighbours to consider for the k-NN
%   - classifier: identification of the chosen classifier
% Return: a fitness array where each value corresponds to the fitness of
% each chromosome of the population + trained ANNs array if the chosen classifier
% is the Multi-Layer Perceptron
function fitness_array = evalPop(mat, pop, pop_size, ...
    nb_features, events, k, classifier)

    for i = (1: pop_size)
        %Build new mat according to features selected
        new_mat = buildSubMat(mat, pop(i, :), nb_features);
            
%         disp('-------------------')
        [accuracy, ~] = classifiers(classifier, new_mat, events, k);
%         accuracy
%         [accuracy, net] = classifiers(classifier, new_mat, events, k, 0);
%         accuracy
%         disp('-------------------')
        fitness_array(i) = accuracy;
%         net_array{i} = net;
    end
end


% Sort the population by fitness (descending order) and keep a record of
% the best fitness achieved for this generation
% Params: 
%   - pop: population of chromosomes
%   - fitness_array: array of fitnesses for one generation
%   - fitness_evolution: array of best fitness per generation
%   - best fitness of the generation
%   - best_features: best chromosome of the population
%   - net_array: array of MLPs for one generation
% Return: 
%   - sorted population
%   - updated array of best fitness per generation
%   - (updated) best fitness
%   - best_features: best individual of the population
%   - best trained ANN per generation (if the MLP is the chosen clasifier)
function [pop, fitness_evolution, best_fitness, best_features, fitness_array] = ...
    sortRecord(pop, fitness_array, fitness_evolution, best_fitness, ...
    best_features)

    [fitness_array, sorted_idx] = sort(fitness_array, 'descend');
    pop = pop(sorted_idx, :);
    if (fitness_array(1) > best_fitness)
        best_features =  pop(1, :);
    end
    best_fitness = max(fitness_array(1), best_fitness);
    fitness_evolution = [fitness_evolution best_fitness];
    
%     pop = [pop fitness_array'];
%     
%     [~, idx] = sort(pop(:, size(pop, 2)), 'descend');
%     
%     pop = pop(idx, :); 
%     
%     if (pop(1, size(pop, 2)) > best_fitness)
%         best_fitness = pop(1, size(pop, 2));
%         best_features =  pop(1, 1:size(pop, 2)-1);
%     end
%     
%     fitness_evolution = [fitness_evolution best_fitness];
%     fitness_array = pop(:, size(pop, 2));
%     pop = pop(:, 1 : size(pop, 2)-1);
end
        

% Population reproduction
% Params: 
%   - pop: population of chromosomes
%   - pop_size: size of the population
%   - crossover_rate
%   - mutation_rate 
%   - fitness_array
%   - nb_features
%   - nb_elites
% Return: new population
function pop = reproduction(pop, pop_size, crossover_rate, mutation_rate, ...
            fitness_array, nb_features, nb_elites)
    new_pop = zeros(pop_size, nb_features);
    %Elitism
    for j = (1: nb_elites)
        new_pop(j, :) = pop(j, :);
    end
    
    % For each individual
    i = nb_elites+1;
    while (i < pop_size)
        % 2.1 Selection: tournament selection
        %[parent1, parent2] = tournamentSelection(pop, pop_size, 4, fitness_array);
        [parent1, parent2] = rouletteWheelSelection(pop, fitness_array, fitness_array(end));
        
        % 2.2 Crossover (one point)
        child1 = parent1;
        child2 = parent2;
        if (rand()<=crossover_rate)
            [child1, child2] = onePointCrossover(parent1, parent2, nb_features);
            %[child1, child2] = twoPointCrossover(parent1, parent2, nb_features);
            %[child1, child2] = uniformCrossover(parent1, parent2, nb_features);
        end

        % 2.3 Mutation
        child1 = mutation(child1, mutation_rate, nb_features);
        child2 = mutation(child2, mutation_rate, nb_features);
        
        new_pop(i, :) = child1;
        new_pop(i+1, :) = child2;
        i = i + 2;
    end
    
    %Update the population
    pop = new_pop;    
end


function [chromosome1, chromosome2] = rouletteWheelSelection(pop, fitnesses, worst_fitness)
    P = fitnesses/worst_fitness;
    P = P/sum(P);
    
    id1 = RouletteWheelSelection(P);
    id2 = id1;
    while (isequal(id1, id2))
        id2 = RouletteWheelSelection(P);
    end
    
    chromosome1 = pop(id1, :);
    chromosome2 = pop(id2, :);
end
   
