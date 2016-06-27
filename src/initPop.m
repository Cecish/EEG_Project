function pop = initPop(pop_size, nb_features )
    
    pop = zeros(pop_size, nb_features);

    for i = (1: pop_size)
        % Creation of random chromosomes
        pop(i, :) = randi([0 1], [1 nb_features]);
    end
end

