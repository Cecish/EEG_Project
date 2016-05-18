function accuracy = OneNearestNeighbour(matrix_min_max, height, width)

    res_accuracy = 0.0;
    
    % For each record
    for i = (1:height)
        ref_target = matrix_min_max(i, width -1);
        
        closest_record = closestNeighbour(i, matrix_min_max, height, width);
        
        %Nearest neighbour has the same class value that the current record
        if isequal(closest_record(width-1), ref_target)
                res_accuracy = res_accuracy + 1;
        end
    end
    
    accuracy = (100*res_accuracy)/height;
end

function closest_neighbour = closestNeighbour(index_record, matrix, height, width)
    new = 0;
    other_distance = 0.0;
    closest_neighbour = [];
    
    for i = (1:width)

        if ~(isequal(i, index_record))

            for j = (1:width-1) %Don't consider the target class!
                other_distance = other_distance + ( (matrix(index_record, j)...
                    - matrix(i, j))^2);
            end

            if isequal(new, 0)
                min_distance = other_distance + 1;
                new = 1;
            end

            %Update best distance && corresponding record
            if (min_distance > other_distance)

                min_distance = other_distance;
                closest_neighbour = matrix(i, :);
            end

            other_distance = 0.0;
        end
    end
end

