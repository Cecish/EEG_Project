function [ hand, k, path_data_file, name_data_file, level, wavelet,...
    classifier, device ] = readParameterFile()

    str_hand = 'HAND: Rehabilitation on the right hand';
    str_device = 'Emotiv EPOC';
    class_name = 'Not supported by this MATLAB program. An error will be caught';

    fileID = fopen('input_file.txt');
    parameters = textscan(fileID,'%d %s %s %d %d %d %s %d','delimiter','\n');
    fclose(fileID);
    
    hand = parameters{1};
    if isequal(hand, 2)
       str_hand = 'HAND: Rehabilitation on the left hand';
    end
    display(str_hand);
    
    path_data_file = parameters{2}{1};  
    display(['Data file path name: ', path_data_file]);
    
    name_data_file = parameters{3}{1};
    display(['Data file name: ', name_data_file]);
    
    classifier = parameters{4};
    switch classifier
        case 1
            class_name = 'k-NN';
        case 2
            class_name = 'SVM';
        case 3
            class_name = 'MLP';
    end
    display(['Classifier name: ', class_name]);
    
    k = parameters{5};
    display(['Number of nearest neighbours: ', num2str(k)]);
    
    level = parameters{6};
    display(['Number of level for the DWT: ', num2str(level)]);
    
    wavelet = parameters{7}{1};
    display(['Mother wavelet for the DWT: ', wavelet]);
    
    device = parameters{8};
    if isequal(device, 2)
       str_device = 'BrainAmp actiCAP';
    end
    display(str_device);
end

