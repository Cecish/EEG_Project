% Menu for letting the user specify his/her parameters
% Return the parameters the user specified: 
%   - hand: hand which wears the robotic glove (either '1'-right or '0'-left)
%   - k: number of nearest neighbours
%   - path_data-file: path of the raw EEG data file
%   - name-data_file: name of the raw EEG data file
%   - level: For the discrete wavelet wavelet transform
%   - wavelet: For the discrete wavelet wavelet transform
%   - classifier
function [ hand, k, path_data_file, name_data_file, level, wavelet,...
    classifier, device ] = menu()
    k = 0;
    
    disp([' ---------------------------------------------------------- ']);
    disp(['|                        Parameters                        |']);
    disp([' ---------------------------------------------------------- ']);
    hand = str2num(input('Rehabilitation on the\n1. right hand\n2. left hand\n','s')); 
    path_data_file = input(strcat('Path of the raw EEG data file?', ... 
        ...%'(H:\\MasterProject\\RawEEGData\\actiCAP\\SAL01\\Online\\)\n'),'s'); 
        '(H:\\MasterProject\\RawEEGData\\EPOC\\No_Glove\\MV\\Online\\)\n'),'s'); 
    name_data_file = input(strcat('Name of the raw EEG data file? ', ...
        ...%' (motor-imagery-csp-4-online-[2015.11.01-11.54.07].vhdr)\n'), 's');
        ' (motor-imagery-csp-4-online-[2016.06.08-18.53.09].vhdr)\n'), 's');
    classifier = str2num(input(strcat('Which classifier?\n 1. kNN\n 2.', ...
        ' SVM\n 3. Multi-Layer Perceptron\n'), 's'));
    % Extra parameter for the k-NN classifier
    if isequal(classifier, 1)
       k = str2num(input('How many nearest neighbours to consider?\n', 's')); 
    end
    level = str2num(input('How many level for the DWT? (recommnded: 5)\n', 's'));
    wavelet = input('Mother wavelet for the DWT? (recommnded: sym1)\n', 's')
end
