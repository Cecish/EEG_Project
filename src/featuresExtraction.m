% Performs the feature extrtaction
% Params: 
%   - mat_X: [records x channels] data matrix
%   - level: decomposition level of the discrete wavelet transform
%   - nb_events: number of events/trials
%   - wavelet
%   - s_rate: sampling rate
%   - signal_length: number of columns for the EEG recorded data (2D matrix)
%   - desired_nb_channels: desired number of optimal channels
% Return: 
%   - mat_features: [40 x nb_channels*level*nb_features] matrix, for each event
%   - ratio_cat_feat: number of features extrcated in each category of
%   features (DWT, AR, PSD, time domain features)
function [mat_features, ratio_cat_feat] = featuresExtraction(mat_X, level, ...
    nb_events, wavelet, s_rate, signal_length, desired_nb_channels)

    ar_order = 6; %Recommended value, as seen in the litterature
    length_block = size(mat_X,1)/40; %number of records for one trial (40 trials in total)
    
    % For each trial/event
    for i = (1: nb_events)
%         temp = ((i-1)*82 +1);
        temp = ((i-1)*length_block +1);
        row = [];

        %For each channel
        for k = (1: desired_nb_channels)
            % 1. #### Discrete wavelet transform ####
            [c,l] = wavedec(mat_X(temp:temp+length_block-1,k),level, wavelet);
%             [c,l] = wavedec(mat_X(temp:temp+81,k),level, wavelet);
            
            %Extract details and approximation coefficients
            [cD, ~] = extractCoeff(level, c, l, wavelet);
            
            % 2. #### Auto-regressive model (AR) ####
            ar_res = ar(mat_X(temp:temp+length_block-1,k), ar_order);
            
            % 3. #### Power Spectrum Analysis ####
            [spectra, freqs] = extractPSDFeatures(s_rate, signal_length, ...
                mat_X(temp:temp+length_block-1, k)); % Focusing on channel k
            psd_features = extractPowerSpectrumFeatures(spectra, freqs);
            
            % 4. #### Time domain features (with Hjorth parameters) ####
            mav = meanabs(mat_X(temp:temp+length_block-1,k)); %Mean absolute value
            mean_val = mean(mat_X(temp:temp+length_block-1,k)); %Mean
            std_val = std(mat_X(temp:temp+length_block-1,k)); %Standard deviation
            
            %Build features row associated to the current trial (exclude cD{1})
            td_features = [mav mean_val std_val];
            row_temp = buildRowFeatures(level, cD, ar_res.A(2:end), ...
                psd_features, td_features);
            
            row = [row row_temp];
        end
        % Update the result matrix with the features extracted for the
        % current trial
        mat_features(i, :) = row;
    end
    ratio_cat_feat = [14 ar_order 10 3]; %[DWT AR PSD Time domain features]
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
%     cD = cell(1, 5); %Initialisation

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
%   - ar_coeffs: AR coefficients
%   - psd_features: Power spectrum analysis features
%   - td_features: Time domain features
% Return: Row features associated to a specific trial
function row = buildRowFeatures(level, cD, ar_coeffs, psd_features, td_features)

    %Initialisation
    shannon_entropy = zeros(1, level-1);
    log_energy = shannon_entropy;
    mean_val1 = shannon_entropy;
    std_val1 = shannon_entropy;
    rms = shannon_entropy;
    mav = shannon_entropy;
    ieeg = shannon_entropy;
    ssi = shannon_entropy;
    var = shannon_entropy;
    aac = shannon_entropy;
    mini = shannon_entropy;
    maxi = shannon_entropy;
    meanVal = shannon_entropy;
    stdVal = shannon_entropy;

    % Features associated to the detail coefficents
    for j = (2: level)
        shannon_entropy(j-1) = wentropy(cD{j}', 'shannon');
        log_energy(j-1) = wentropy(cD{j}', 'log energy');
        [ey, ~] = energyop(cD{j}', 0); %No plot
        mean_val1(j-1) = mean(ey);
        std_val1(j-1) = std(ey);
        rms(j-1) = rootMeanSquare(cD{j}');
        mav(j-1) = meanAbsoluteValue(cD{j}');
        ieeg(j-1) = integratedEEG(cD{j}');
        ssi(j-1) = simpleSquareIntegral(cD{j}');
        var(j-1) = varianceEEG(cD{j}');
        aac(j-1) = averageAmplitudeChange(cD{j}');
        mini(j-1) = min(cD{j});
        maxi(j-1) = max(cD{j});
        meanVal(j-1) = mean(cD{j});
        stdVal(j-1) = std(cD{j});
    end
    
    row = [rms mav ieeg ssi var aac mini maxi meanVal stdVal shannon_entropy...
        log_energy mean_val1 std_val1 ar_coeffs psd_features td_features];
%     row = [rms mav ieeg ssi var aac mini maxi meanVal stdVal ar_coeffs ...
%         psd_features td_features];
end


% Power Spectrum analysis: extraction of features
% Params: 
%   - spectra: power spectrum of a specific channel
%   - freqs: 
% Return: array of features extracted from the power spectrum analysis
function psd_features = extractPowerSpectrumFeatures(spectra, freqs)
        
    %Extraction of features for the five frequency bands
    %[2 - 4] Hz
    band1_idx = find(freqs>=2 & freqs<=4);  
    band1_power = mean(spectra(band1_idx)); %Mean
    band1_std = std(spectra(band1_idx)); %Std
    %Standard deviation
    %[4 - 8] Hz
    band2_idx = find(freqs>=4 & freqs<=8);  
    band2_power = mean(spectra(band2_idx)); %Mean
    band2_std = std(spectra(band2_idx)); %Std
    %[8 - 12] Hz
    band3_idx = find(freqs>=8 & freqs<=12); 
    band3_power = mean(spectra(band3_idx)); %Mean
    band3_std = std(spectra(band3_idx)); %Std
    %[12 - 18] Hz
    band4_idx = find(freqs>=12 & freqs<=18); 
    band4_power = mean(spectra(band4_idx)); %Mean
    band4_std = std(spectra(band4_idx)); %Std
    %[18 - 30] Hz
    band5_idx = find(freqs>=18 & freqs<=30);
    band5_power = mean(spectra(band5_idx)); %Mean
    band5_std = std(spectra(band5_idx)); %Std
       
    %Update res = matrix where each line gathers these extracted
    %features for a specific channel
    psd_features = [band1_power band1_std band2_power band2_std band3_power....
        band3_std band4_power band4_std band5_power band5_std];
end


% Params:
%   - s_rate: sampling frequency
%   - signal_length: length of signal
%   - x: acquisition data from a specific channel
% Return: 
%   - spectra: power spectrum of a specific channel
%   - freqs: frequences
function [spectra,freqs] = extractPSDFeatures(s_rate, signal_length, x)

    NFFT = 2^nextpow2(signal_length);     % Next power of 2 from length of x
    NOVERLAP = 0;
    WINDOW = 512;

    %Matlab pwelch function
    [spectra,freqs] = pwelch(x,WINDOW,NOVERLAP,NFFT,s_rate);
%    [spectra,freqs] = spectopo(mat_X(temp:temp+length_block-1,k),...
%         0, eeg.srate, 'nfft', 1024); %
            
end
