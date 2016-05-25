% Menu for letting the user specify his/her parameters
% Return the parameters the user specified: 
%   - hand: hand which wears the robotic glove (either 'r'-right or 'l'-left)
%   - k: number of nearest neighbours
%   - path_data-file: path of the raw EEG data file
%   - name-data_file: name of the raw EEG data file
%   - level: For the discrete wavelet wavelet transform
%   - wavelet: For the discrete wavelet wavelet transform
%   - classifier
function [ hand, k, path_data_file, name_data_file, level, wavelet,...
    classifier ] = menu()
    k = 0;
    
    disp([' ---------------------------------------------------------- ']);
    disp(['|                        Parameters                        |']);
    disp([' ---------------------------------------------------------- ']);
    hand = input('Rehabilitation on the right(r) or left(l) hand?\n','s'); 
    path_data_file = input(strcat('Path of the raw EEG data file?', ... 
        '(H:\\MasterProject\\RawEEGData\\EPOC\\acquisition\\ToUse-Alistair\\)\n'),'s'); 
    name_data_file = input(strcat('Name of the raw EEG data file? ', ...
        ' (motor-imagery-csp-1-acquisition-[2016.03.18-16.55.10]_alistair.vdhr)\n'), 's');
    classifier = str2num(input('Which classifier?\n 1: kNN\n 2: SVM\n', 's'));
    if isequal(classifier, 1)
       k = str2num(input('How many nearest neighbours to consider?\n', 's')); 
    end
    level = str2num(input('How many level for the DWT? (recommnded: 5)\n', 's'));
    wavelet = input('Mother wavelet for the DWT? (recommnded: sym1)\n', 's')
end
