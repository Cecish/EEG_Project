% Functions are fields of a struct (holder)
% Return function handles to local functions
function funs = featureSelection
    funs.featSelection = @featSelection; %Features selection
    funs.chanReduction = @chanReduction; %Best channels subset selection
end

% Feature selection
% Params: 
%   - mat: initial matrix with all the features previously extracted
%   - nb_features: total number of features
%   - k: number of nearest neighbours to consider with the k-NN classifier
%   - events: known outcome associated to each trial
%   - classifier: id of the chosen classifier 
%   - do_plot: plotting or not accuracy/fitness/penalty evolutions (0: No, 1: Yes)
% Return: 
%   - res: sub matrix (with a reduced number of features) that led to the
%     best accuracy score
%   - best_features: subset of features that led to the best accuracy score
function [res, best_features] = featSelection( mat, nb_features, k,...
    events, classifier, do_plot)
    
    % GA parameters
    pop_size = 20;
    max_generations = 50;
    nb_elites = 2;
    crossover_rate = 0.8;
    mutation_rate = 0.1;
    
    % Parameters initialisation
    best_features = [];
    it = 0;
    best_model = 0;
    best_fitness = -Inf;
    fitness_evolution = [];
    accuracy_evolution = fitness_evolution;
    penalty_evolution = fitness_evolution;
    
    % Population initialisation
    [pop, fitness_array, accuracy_array, penalty_array, fitness_evolution, ...
        best_fitness, best_features, accuracy_evolution, penalty_evolution, net_array, best_model] = ...
        popInit(pop_size, nb_features, fitness_evolution, best_fitness, ...
        best_features, accuracy_evolution, penalty_evolution, mat, events,...
        k, classifier, best_model);
        
    %Loop until the termination condition is reached
    while ((it < max_generations) && (best_fitness < 100.00))
        disp(['#### EPOCH n° ', num2str(it)]);

        % ####1: Create a new population
        [pop, fitness_array, accuracy_array, penalty_array, net_array] = reproduction(...
            pop, pop_size, crossover_rate, mutation_rate, fitness_array, ...
            nb_features, nb_elites, accuracy_array, penalty_array, mat, ...
            classifier, events, k, net_array);
        
        %####2: Sort by population by fitness value (descending order)
        [pop, fitness_evolution, best_fitness, best_features, best_model, fitness_array,...
            accuracy_evolution, penalty_evolution] = sortRecord(pop, ...
            fitness_array, fitness_evolution, best_fitness, best_features,...
            accuracy_evolution, penalty_evolution, nb_features, net_array, best_model);
    
        %####3: Truncate the population (only pop_size individuals per population)
        pop = pop(1:pop_size, :);
        fitness_array = fitness_array(1:pop_size);
        accuracy_array = accuracy_array(1:pop_size);
        penalty_array = penalty_array(1:pop_size);
        
        it = it + 1; % Increment
        display(['Iteration: ',num2str(it),': Best fitness = ', ...
            num2str(fitness_evolution(it+1)), ', Best accuracy = ', ...
            num2str(accuracy_evolution(it+1)), ', Best penalty = ', ...
            num2str(penalty_evolution(it+1))]);
    end

    %Save the fitness evolution in a file
%     dlmwrite('fitness_evol.csv', fitness_evolution', 'delimiter', ';');
    
    disp(['Best features: ', mat2str(best_features)]);
    
    %Build matrix with the best features
    res = buildSubMat(mat, best_features, nb_features);
    
    if (do_plot)
        % Plotting the evolution of the best fitness value, the "best "accuracy
        % score (here, best in terms of the corresponding fitness value) and the
        % "best" penalty term (after penalty and being part of the best fitness
        % value) over the generations 
        abscissa = 1:it;
        
        figure
        evolution_graph = plot(abscissa, fitness_evolution(2:end), abscissa, ...
            accuracy_evolution(2:end), abscissa, penalty_evolution(2:end));
        grid on;
        title('Best fitness/accuracy/penalty evolution over the generations')
        xlabel('Generations')
        ylabel('Fitness/Accuracy/Rest after penalty')
        legend(evolution_graph, {'Fitness evolution', ...
        'Evolution of the accuracy associated with the best fitness value',...
        'Evolution of the best remainder after penalties, associated with the best fitness value'}, ...
        'Location','Best','FontSize',8);
    end
    
end


% Randomly create a chromosome
% Param: nb_features: number of genes for the chromosome
% Return: new chromosome
function chromosome = createIndividual(nb_features)
    chromosome = zeros(1, nb_features); %Initialisation

    for i = (1: nb_features)
        if (rand()<0.5)
            chromosome(i) = 0;
        else
            chromosome(i) = 1;
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
    
    % Count number of selected features (number of 1 in the binary array)
    nb_selected_feat = sum(individual == 1); %Number of features selected
    new_mat = zeros(size(mat, 1), nb_selected_feat); %Initialisation
    
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
%   - fitness_array: array of features
% Return: the two parents (chromosomes) selected
function [parent1, parent2] = tournamentSelection(pop, pop_size, nb, fitness_array)
    
    temp = zeros(nb, size(pop, 2)+1); %Initalisation

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
    
    % Children initialisation
    child1 = zeros(1, length);
    child2 = child1;

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
    
    %Children initialisation
    child1 = zeros(1, length);
    child2 = child1;

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
%   - length: length of the chromosomes
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
function pop = createPop(pop_size, nb_features)

    pop = zeros(pop_size, nb_features); %Initialisation

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
% Return: 
%   - fitness_array: a fitness array where each value corresponds to the 
% fitness of each chromosome of the population (accuracy + penalty term)
%   - accuracy_array: array of the accuracies gotten by each individual of 
% the population
%   - penalty_array: array of penalties. Its values are the number of
%   fetures not selected for a specific solution. It is thus
%   intended/encouraged that the values will increase over the generations
function [fitness_array, accuracy_array, penalty_array, net_array] = evalPop(mat, pop,...
    pop_size, nb_features, events, k, classifier)

    %Arrays initialisation
    fitness_array = zeros(1, pop_size);
    accuracy_array = fitness_array;
    penalty_array = fitness_array;

    for i = (1: pop_size)
        %Build new mat according to features selected
        new_mat = buildSubMat(mat, pop(i, :), nb_features);
            
        % Calculation of the accuracy
        [accuracy, net_array{i}] = classifiers(classifier, new_mat, events, k);
        accuracy_array(i) = accuracy;
        % Calculation of the fitness_value
        % AND keep note of the penalty value of each solution of the population
        nb_extracted_feat = sum(pop(i, :) == 1); %Number of features selected
        penalty_array(i) = (nb_features - nb_extracted_feat)/nb_features*100;
        fitness_array(i) = 0.5*accuracy + 0.5*penalty_array(i);
    end
end


% Sort the population by fitness (descending order) and keep a record of
% the best fitness achieved for this generation
% Params: 
%   - pop: population of chromosomes
%   - fitness_array: array of fitnesses for one generation
%   - fitness_evolution: array of best fitness per generation
%   - best_fitness: best fitness of the generation
%   - best_features: best chromosome of the population
%   - accuracy_evolution: evolution of the "best" accuracies over the
%   generations (=accuracy of the best solution)
%   - penalty_evolution: penalty_evolution: evolution of the "best" remainders after penalty 
%   over the generations
%   - nb_features: number of genes for each individual 
% Return: 
%   - pop: sorted population
%   - fitness_evolution: updated array of best fitness per generation
%   - best_fitness: (updated) best fitness
%   - best_features: best individual of the population
%   - fitness_array: updated array of fitnesses for one generation
%   - accuracy_evolution: evolution of the "best" accuracies over the
%   generations (=accuracy of the best solution) (updated)
%   - penalty_evolution: penalty_evolution: evolution of the "best" remainders after penalty 
%   over the generations (updated)
function [pop, fitness_evolution, best_fitness, best_features, best_model, ...
    fitness_array, accuracy_evolution, penalty_evolution] = sortRecord(...
    pop, fitness_array, fitness_evolution, best_fitness,...
    best_features, accuracy_evolution, penalty_evolution, nb_features, net_array, best_model)

    [fitness_array, sorted_idx] = sort(fitness_array, 'descend');
    pop = pop(sorted_idx, :);
 
    if (fitness_array(1) > best_fitness)
        best_features =  pop(1, :);
% save('truc', 'net_array');
        best_model = net_array(sorted_idx(1));
    end
    best_fitness = max(fitness_array(1), best_fitness);
    fitness_evolution = [fitness_evolution best_fitness];
    nb_extracted_feat = sum(best_features == 1);
%     best_accuracy = 2*best_fitness - ((nb_features - nb_extracted_feat)/nb_features);
    best_accuracy = (best_fitness - 0.5*((nb_features - nb_extracted_feat)/nb_features)*100)/0.5;
    best_penalty = (best_fitness - 0.5*best_accuracy)/0.5; %Smallest penalty

    accuracy_evolution = [accuracy_evolution best_accuracy];
    penalty_evolution = [penalty_evolution best_penalty];
end
        

% Population reproduction
% Params: 
%   - pop: population of chromosomes
%   - pop_size: size of the population
%   - crossover_rate
%   - mutation_rate 
%   - fitness_array: array of fitnesses for one generation
%   - nb_features
%   - nb_elites
%   - accuracy_array: array of accuracies for one generation
%   - penalty_array: array of "penalties" for one generation
%   - mat: matrix of features
%   - classifier: identification of the chosen classifier
%   - events: events associated to each row of mat
%   - k: number of neighbours to consider for the k-NN
% Return: 
%   - pop: new population
%   - fitness_array_bis: updated array of fitnesses corresponding to 40
%   individuals (20 fitnesses associated to the original population and 20
%   fitnesses associated to the new population of offsprings concatenated
%   to the old population for now)
%   - accuracy_array: updated array of accuracies for one generation
%   - penalty_array: updated array of "penalties" for one generation
function [pop, fitness_array_bis, accuracy_array, penalty_array, net_array] = ...
    reproduction(pop, pop_size, crossover_rate, mutation_rate, fitness_array,...
    nb_features, nb_elites, accuracy_array, penalty_array, mat, classifier,...
    events, k, net_array)

    new_pop = zeros(pop_size, nb_features); %Initialisation
    
    fitness_array_bis = fitness_array;
    
    %Elitism
    new_pop(1:nb_elites, :) = pop(1:nb_elites, :);
    
    % For each individual
    i = nb_elites+1;
    while (i < pop_size)
        % 2.1 Selection: roulette wheel selection (comment/uncomment)
%         [parent1, parent2] = tournamentSelection(pop, pop_size, 4, fitness_array);
        [parent1, parent2] = rouletteWheelSelection(pop, fitness_array, fitness_array(end));
        
        % 2.2 Crossover, one point (commet/uncomment)
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
        
        % Evaluation of the children. Doing it there avoid evaluating the
        % whole population later (and re-evaluate already existing and thus
        % evaluated individuals)
        accuracy_array(pop_size+1:pop_size+nb_elites) = accuracy_array(1:nb_elites);
        penalty_array(pop_size+1:pop_size+nb_elites) = penalty_array(1:nb_elites);
        fitness_array_bis(pop_size+1:pop_size+nb_elites) = fitness_array(1:nb_elites);
        net_array(pop_size+1:pop_size+nb_elites) = net_array(1:nb_elites);
        
        % Processing of child n°1
        [fitness_array_bis, accuracy_array, penalty_array, net_array] = evalChild(...
            accuracy_array, penalty_array, fitness_array_bis, mat, nb_features,...
            child1, classifier, events, k, pop_size, i, net_array);
        
        % Processing of child n°2
        [fitness_array_bis, accuracy_array, penalty_array, net_array] = evalChild(...
            accuracy_array, penalty_array, fitness_array_bis, mat, nb_features,...
            child2, classifier, events, k, pop_size, i+1, net_array);
        
        i = i + 2;
    end
    
    %Update the population 
    pop = [pop; new_pop];
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
   
% Best channels subset selection
% Params: 
%   - mat: initial matrix NxM where N is the number of samples and M is the
%   number of original channels
%   - desired_nb_channels: number of best channels to keep for the analysis
% Return: 
%   - reduced_mat: The reduced matrix Nxdesired_nb_channels
%   - best_channels: the first desired_nb_channels optimal channels
function [reduced_mat, best_channels, worst] = chanReduction(mat, ex_events_Y, desired_nb_channels)

    selectedIndices = mrmr_mid_d(mat, ex_events_Y, size(mat, 2))
   
    best_channels = selectedIndices(1:desired_nb_channels);
    worst_channels = selectedIndices(desired_nb_channels+1:end);
    
    reduced_mat = mat(:, best_channels);
    worst = mat(:, worst_channels);
end


% Population initialisation. This includes a random instanciation of each
% individual of the population, their evaluation with a specific classifier
% and this population's sorting by fitness value (composed of an accuracy 
% and a penalty term)
% Params: 
%   - pop_size: size of the population to create
%   - nb_features: number of genes for each individual
%   - fitness_evolution: evolution of best fitnesses over the generations
%   - best_fitness: best fitness of the generation
%   - best_features: best chromosome of the population
%   - accuracy_evolution: evolution of the "best" accuracies over the
%   generations (=accuracy of the best solution)
%   - penalty_evolution: evolution of the "best" remainders after penalty 
%   over the generations
%   - mat: matrix of features
%   - events: events associated to each row of mat
%   - k: number of neighbours to consider for the k-NN
%   - classifier: identification of the chosen classifier
% Return: 
%   - pop: new population (instanciated, evaluated and sorted)
%   - fitness_array: array of fitnesses for one generation
%   - accuracy_array: array of accuracies for one generation
%   - penalty_array: array of "penalties" for one generation
%   - fitness_evolution: array of best fitness per generation
%   - best_fitness: best fitness of the generation (updated)
%   - best_features: best chromosome of the population (updated)
%   - accuracy_evolution: evolution of the "best" accuracies over the
%   generations (=accuracy of the best solution) (updated)
%   - penalty_evolution: evolution of the "best" remainders after penalty 
%   over the generations (updated)
function [pop, fitness_array, accuracy_array, penalty_array, fitness_evolution, ...
        best_fitness, best_features, accuracy_evolution, penalty_evolution, net_array, best_model] = ...
        popInit(pop_size, nb_features, fitness_evolution, best_fitness, ...
        best_features, accuracy_evolution, penalty_evolution, mat, events,...
        k, classifier, best_model)
    
    % Population creation
    pop = createPop(pop_size, nb_features);
    
    % Initial population evaluation
    [fitness_array, accuracy_array, penalty_array, net_array] = evalPop(mat, pop,...
        size(pop, 1), nb_features, events, k, classifier);
    
    % Sort the population by fitness and record the best fitness of the generation
    [pop, fitness_evolution, best_fitness, best_features, best_model, fitness_array,...
        accuracy_evolution, penalty_evolution] = sortRecord(pop, fitness_array, fitness_evolution,...
        best_fitness, best_features, accuracy_evolution, penalty_evolution, nb_features, net_array, best_model);
    
end


% Evaluation of the children obtained after crossover + mutation
% Params: 
%   - accuracy_array: array of accuracies for one generation
%   - penalty_array: array of "penalties" foor one generation
%   - fitness_array_bis: array of fitnesses for one generation
%   - mat: matrix of features
%   - nb_features: number of genes for each individual
%   - child: individual to evaluate
%   - classifier: identification of the chosen classifier
%   - events: events associated to each row of mat
%   - k: number of neighbours to consider for the k-NN
%   - pop_size: size of a population
%   - i: index of the new individual currently being processed
% Return: 
%   - fitness_array_bis: array of fitnesses for one generation (updated)
%   - accuracy_array: array of accuracies for one generation (updated)
%   - penalty_array: array of "penalties" foor one generation (updated)
function [fitness_array_bis, accuracy_array, penalty_array, net_array] = evalChild(...
            accuracy_array, penalty_array, fitness_array_bis, mat, nb_features,...
            child, classifier, events, k, pop_size, i, net_array)
        
    %Build new mat according to features selected
    new_mat = buildSubMat(mat, child, nb_features);
    % Calculation of the accuracy
    [accuracy, net_array{pop_size+i}] = classifiers(classifier, new_mat, events, k);
    accuracy_array(pop_size+i) = accuracy;
    % Calculation of the penalty value
    nb_extracted_feat = sum(child == 1); %Number of features selected
    penalty_array(pop_size+i) = (nb_features - nb_extracted_feat)/nb_features*100;
    % Calculation of the fitness value
    fitness_array_bis(pop_size+i) = 0.5*accuracy + 0.5*penalty_array(pop_size+i);
        
end