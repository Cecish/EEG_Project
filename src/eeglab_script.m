% Instead of using the EEGLAB GUI
% Params:
%   - path: path of the raw EEG data file
%   - name: name of the raw EEG data file
%   - highpass: highpass filter
%   - notch: notch filter
%   - device: device identifier (1: EPOC, 2: actiCAP)
% Return: a alleeg structure similar to the ALLEEG variable produced by
% eeglab toolbox when using the GUI
function alleeg = eeglab_func(path, name, highpass, notch, device)
    eegh;
    
    % #### 0: Load the data file into EEGLAB toolbox
    EEG = pop_loadbv(path, name);
    EEG.setname = 'raw';
        
    % #### 1: Edit channels location
    EEG = pop_chanedit(EEG, 'lookup');
    alleeg(1) = EEG;
    
    % #### 2: Preprocessing
    % 2.1 Re-reference (average) ---------> 23.8934%
    EEG = pop_reref( EEG, []);
    EEG.setname = 'Alistair re-referencing average';
    alleeg(2) = EEG;
    
    % 2.2 High pass filter (1Hz) ---------> 20.9124%
    %EEG = pop_eegfiltnew(EEG, [], highpass, 424, true, [], 0);
    EEG = pop_eegfiltnew(EEG, [], highpass, 6600, true, [], 0);

    EEG.setname = 'Alistair highpass filter';
    alleeg(3) = EEG;
    
    % 2.3 Notch filter ------------> 20.8672%
    %EEG = pop_eegfiltnew(EEG, notch(1), notch(2), 212, 1, [], 0);
    EEG = pop_eegfiltnew(EEG, notch(1), notch(2), 6600, 1, [], 0);
    EEG.setname = 'Alistair notch filter';
    alleeg(4) = EEG;
    
    %EEG = pop_autobsseog( EEG, [98], [98], 'bsscca', {'eigratio', [1000000]}, 'eog_fd', {'range',[1  4]});
    EEG = pop_autobsseog( EEG, [128], [128], 'sobi', {'eigratio', [1000000]}, 'eog_fd', {'range',[1  5]});
    alleeg(5) = EEG;
    
    % #### X: Epoch extraction
    epochs = extractEpochInterval(device);
    EEG = pop_epoch( EEG, {'S769' 'S770'}, epochs, 'newname', ...
        'Alistair epochs', 'epochinfo', 'yes');
    alleeg(6) = EEG;
    
    %Reject threshold [-100, 100 uV]
    %EEG = pop_eegthresh(EEG,1,[1:14] ,-100,100,0,0.63281,0,0);
    %EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
    %EEG = pop_rejepoch( EEG, [4 5 12 13 14 17 22 23 26 27 28 31 32 33 35 37 40] ,0);
    %EEG.setname='Alistair reject threshold';
    %alleeg(7) = EEG;
end


% Extract the right window interval for epoch extraction, depending on the
% EEG device used: [starting the event 769 or 770 - end of the event]ms: 
% - [0 0.640] for the Emotiv EPOC headset
% - [0 1] for the actiCAP cap
% Parameters: device (1: EPOC, 2: actiCAP)
% Return: the right window interval
function epoch_interval = extractEpochInterval(device)
    
    if isequal(device, 1) %EPOC
        epoch_interval = [0 0.640];
    else %actiCAP
        epoch_interval = [0 1]
    end
end


