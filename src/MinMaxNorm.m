function res_matrix = MinMaxNorm(dataset, height, width)

    minmax_arrays = findMinMax(dataset, height, width);

    for j = (1: width)
        for i = (1 : height)

            if ~isequal(j, width-1)
            	res_matrix(i, j) = (dataset(i, j)-minmax_arrays(1, j))/...
                        (minmax_arrays(2, j)-minmax_arrays(1, j));
            else %Classfield
            	res_matrix(i, j) = dataset(i, j);
            end
        end
    end
end


function minMax_arrays = findMinMax(dataset, height, width)
    
    %Min
    for j = (1: width)
        min = Inf;
        for i = (1:height)
            if (dataset(i, j) < min)
                min = dataset(i, j);
            end
        end
        
        disp(['Minimum of column', num2str(j), ' = ', num2str(min)]);
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
        
        disp(['Maximum of column', num2str(i), ' = ', num2str(max)]);
        minMax_arrays(2, i) = max;
    end
end

