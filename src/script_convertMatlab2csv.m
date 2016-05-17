% Author: Cécile Riquart
% Date 16/05/2016, 14:54

% Converts eeglab data into csv file (whose format is similar to the one
% obtained using the ov2csv scenario)
function convertMatlab2csv(struct_eeg)

    format long; %To be sure that no information is lost during the conversion

    %If EEG.data is a 2D matrix ()
    if isequal(length(size(struct_eeg.data)), 2)
       
        %Build 1st line
        first_line = build_1st_line(struct_eeg);
        
        %Build times + corresponding measures (without sampling rate)
        last_col = zeros(struct_eeg.pnts,1);
        matrix = [struct_eeg.times' struct_eeg.data' last_col];
        matrix(1,struct_eeg.nbchan + 2) = struct_eeg.srate;
        
        %new_l = [char(matrix), last_col];
        %m = cell2mat(last_col);
        %display(last_col);
        %new_mat = [matrix last_col];
        dlmwrite('temp.csv', matrix, 'delimiter', ';');
        %display(size(m));
        %display(size(matrix));
        %C = vertcat(first_line, new_matrix);
        %display(new_mat);
        
        %dlmwrite('temp.csv',new_mat,'delimiter',';');
        
        %Add sampling rate at the end of the second row
        

    end
    
    %If EEG.data is a 3D matrix ()
    
end

%Build 1st line: labels (time + channels names + sampling rate)
function line = build_1st_line(struct_eeg)
    line = 'Time (s);'
    temp = char(struct_eeg.chanlocs.labels);
    temp2 = cellstr(temp);

    for i = 1: struct_eeg.nbchan
        line = strcat(line, temp2(i), ';');
    end
    
    line = strcat(line, 'Sampling Rate');
end