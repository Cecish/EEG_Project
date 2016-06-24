function accuracy = fitnessFunction( original_mat, pop, pop_size, events, classifier, k)

    for i = (1: pop_size)
        new_mat = [];
        new_mat = buildSubMat(original_mat, pop(i, :), size(original_mat, 2));
        accuracy(i) = 100 - classifiers(classifier, new_mat, events, k);
    end
end


function new_mat = buildSubMat(mat, individual, length_id)
    count = 1;

    for i = (1 : length_id)
        if isequal(individual(i), 1) %Features taken for new_mat
            new_mat(:, count) = mat(:, i);
            count = count + 1;
        end
    end
end