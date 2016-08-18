% Instead of using the EEGLAB GUI
% Params:
%   - path: path of the raw EEG data file
%   - name: name of the raw EEG data file
%   - highpass: highpass filter
%   - lowpass: lowpass filter
%   - notch: notch filter
% Return: a alleeg structure similar to the ALLEEG variable produced by
% eeglab toolbox when using the GUI
function alleeg = eeglab_script(path, name, highpass, lowpass, notch, str_device)

    % #### 0: Load the data file into EEGLAB toolbox
    EEG = pop_loadbv(path, name);
    EEG.setname = 'Raw';
        
    % #### 1: Edit channels location
    EEG = pop_chanedit(EEG, 'lookup','H:\MasterProject\EEGLab\eeglab_current\eeglab13_5_4b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp');
    EEG = eeg_checkset( EEG );
%     EEG = pop_chanedit(EEG, 'lookup');
    alleeg(1) = EEG;
    
    % #### 2: Preprocessing    
    % 2.2 High pass filter (1Hz)          
    %EEG = pop_eegfilt( EEG, highpass, 0, [], [0]);
    EEG = pop_eegfiltnew(EEG, highpass, []);
    EEG.setname = 'Highpass filter';
    alleeg(2) = EEG;
    
    %EEG = pop_eegfilt( EEG, 0, lowpass, [], [0], 0, 0, 'fir1', 0);
    EEG = pop_eegfiltnew(EEG, [], lowpass);
    EEG.setname='Lowpass filter';
    alleeg(3) = EEG;
    
    % 2.3 Notch filter
    if isequal(str_device, 2) %actiCAP
        %60Hz power line in Brasil (the actiCAP data has been acquired in Brasil)
        notch(1) = 58;
        notch(2) = 62; 
    end
    EEG = pop_eegfiltnew(EEG, notch(1), notch(2), [], 1, [], []);
    EEG.setname = 'Notch filter';
    alleeg(4) = EEG;
     
    % 2.1 Re-reference (average)
    EEG = pop_reref( EEG, []);
    EEG.setname = 'Re-referencing average';
    alleeg(5) = EEG;
     
    % EOG removal using BSS
    %EEG = pop_autobsseog( EEG, [98], [98], 'bsscca', {'eigratio', [1000000]}, 'eog_fd', {'range',[1  4]});
    if isequal(str_device, 1)
       EEG = pop_autobsseog( EEG, [10], [10], 'sobi', {'eigratio', [1000000]}, 'eog_fd', {'range',[0  5]});
    else
       EEG = pop_autobsseog( EEG, [128], [128], 'sobi', {'eigratio', [1000000]}, 'eog_fd', {'range',[0  5]});
    end
    alleeg(6) = EEG;
    
    % #### 3: Epoch extraction
    %epochs = extractEpochInterval(device);
    epochs = [-1 5];
    EEG = pop_epoch( EEG, {'S769' 'S770'}, epochs, 'newname', ...
        'Epochs extraction', 'epochinfo', 'yes');
    alleeg(7) = EEG;
    
end


