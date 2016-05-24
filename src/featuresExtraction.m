% Performs the feature extrtaction
% Params: 
%   - mat_X: [records x channels] data matrix
%   - level: decomposition level of the discrete wavelet transform
%   - nb_events: number of events/trials
%   - wavelet
% Return: [40 x nb_channels*level*nb_features] matrix, for each event
function mat_features = featuresExtraction(mat_X, level, nb_events, wavelet)
    
    % For each trial/event
    for i = (1: nb_events)
        temp = ((i-1)*82 +1);
        row = [];

        %For each channel
        for k = (1: 14)
            % Discrete wavelet transform
            [c,l] = wavedec(mat_X(temp:temp+81,k),level, wavelet);
            
            %Extract details and approximation coefficients
            [cD, cA] = extractCoeff(level, c, l, wavelet);
            
            %Build features row associated to the current trial
            row_temp = buildRowFeatures(level, cD, cA);

            row = [row row_temp];
        end
        % Update the result matrix with the features extracted for the
        % current trial
        mat_features(i, :) = row;
    end
end


% Discrete wavelet transform
% Params:
%   - mat_X: records x channels data
%   - wavelet_level: 
%   - lowpass: low pass filter
%   - highpass: high pass filter
% Return: 
%   - A: approximation coefficients
%   - D: detail coefficients
%function [cA cD] = dwt_func(mat_X, wavelet_level, lowpass, highpass, name)
%
%    % For eeg channel data
%    %for i = (1: size(mat_X, 2))
%        % Get the coefficients for level 1
%        %[cA,cD{1}] = dwt(mat_X(:, 1)', lowpass(1), highpass(1));
%        [cA,cD{1}] = dwt(mat_X(:, 1)', name);
%
%        % Continue to cascade to get additional coefficients at each level
%        for n = 2:wavelet_level
%            %[cA,cD{n}] = dwt(cA, lowpass(n), highpass(n));
%            [cA,cD{n}] = dwt(cA, name);
%        end
%    %end
%end


% Root Mean Square
% Params: data array
% Return the corresponding root mean square
function rms = rootMeanSquare(data_array)
    temp = 0.0;

    for i = (1: length(data_array))
        temp = temp + data_array(1).^2;
    end
    
    rms = sqrt((1/length(data_array)) * temp);
end


% Mean Absolute Value
% Params: data array
% Return the corresponding mean absolute value
function mav = meanAbsoluteValue(data_array)
    temp = 0.0;
    
    for i = (1: length(data_array))
        temp = temp + abs(data_array(i));
    end
    
    mav = (1/length(data_array)) * temp;
end


% Integrated EEG
% Params: data array
% Return the corresponding integrated EEG
function ieeg = integratedEEG(data_array)
    ieeg = 0.0;
    
    for i = (1:length(data_array))
        ieeg = ieeg + abs(data_array(i));
    end
end


% Simple Square Integral
% Params: data array
% Return the corresponding simple square integral
function ssi = simpleSquareIntegral(data_array)
    ssi = 0.0;
    
    for i = (1:length(data_array))
       ssi = ssi + (abs(data_array(i)).^2); 
    end
end


% Variance of EEG
% Params: data array
% Return the corresponding variance of EEG
function var = varianceEEG(data_array)
    temp = 0.0;

    for i = (1:length(data_array))
        temp = temp + data_array(i).^2;
    end
    
    var = (1/(length(data_array)-1)) * temp;
end


% Average Amplitude Change
% Params: data array
% Return the corresponding average amplitude change
function aac = averageAmplitudeChange(data_array)
    temp = 0.0;
    
    for i = (1: length(data_array)-1)
        temp = temp + abs( data_array(i+1) - data_array(i) );
    end
    
    aac = (1/length(data_array)) * temp;
end


% Extraction of the detail and approximation coefficients from the wavelet
% decomposition level and bookkeeping vectors which resulted from the DWT
% Params: 
%   - level: decomposition level of the discrete wavelet transform
%   - c: wavelet decomposition vector
%   - l: bookkeeping vector
%   - wavelet used for the discrete wavelet transfotm
% Return: 
%   - cD: Detail coeffients structure
%   - cA: Approximation coefficients
function [cD, cA] = extractCoeff(level, c, l, wavelet)
    %Detail coefficients
    for i = (1:level)
       cD{i} = detcoef(c,l, i);
    end
    
    %Approximation coefficients
    cA = appcoef(c,l, wavelet, level);
end


% Build a row of features associated to a specific trial
% Params:
%   - level: decomposition level of the discrete wavelet transform
%   - cD: Detail coeffients structure
%   - cA: Approximation coefficients
% Return: Row features associated to a specific trial
function row = buildRowFeatures(level, cD, cA)
    % Features associated to the detail coefficents
    for j = (1: level)
        rms(j) = rootMeanSquare(cD{j}');
        mav(j) = meanAbsoluteValue(cD{j}');
        ieeg(j) = integratedEEG(cD{j}');
        ssi(j) = simpleSquareIntegral(cD{j}');
        var(j) = varianceEEG(cD{j}');
        aac(j) = averageAmplitudeChange(cD{j}');
    end
    
    %Features associated to the approximation coeffients
    rms(level) = rootMeanSquare(a4);
    mav(level) = meanAbsoluteValue(a4);
    ieeg(level) = integratedEEG(a4);
    ssi(level) = simpleSquareIntegral(a4);
    var(level) = varianceEEG(a4);
    aac(level) = averageAmplitudeChange(a4);
    
    row = [rms mav ieeg ssi var aac];
end
