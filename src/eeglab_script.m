% Instead of using the EEGLAB GUI
% Params:
%   - path: path of the raw EEG data file
%   - name: name of the raw EEG data file
%   - highpass: highpass filter
%   - notch: notch filter
% Return: a alleeg structure similar to the ALLEEG variable produced by
% eeglab toolbox when using the GUI
function alleeg = eeglab_script(path, name, highpass, notch)
    
    % #### 0: Load the data file into EEGLAB toolbox
    EEG = pop_loadbv(path, name);
    EEG.setname = 'Raw';
        
    % #### 1: Edit channels location
    EEG = pop_chanedit(EEG, 'lookup');
    alleeg(1) = EEG;
    
    % #### 2: Preprocessing    
    % 2.2 High pass filter (1Hz)          
    EEG = pop_eegfilt( EEG, highpass, 0, [], [0]);
    EEG.setname = 'Highpass filter';
    alleeg(2) = EEG;
     
    % 2.3 Notch filter
    EEG = pop_iirfilt(EEG,notch(1),notch(2),0,1);
    EEG.setname = 'Notch filter';
    alleeg(3) = EEG;
     
    % 2.1 Re-reference (average)
    EEG = pop_reref( EEG, []);
    EEG.setname = 'Re-referencing average';
    alleeg(4) = EEG;
     
    % EOG removal using BSS
    %EEG = pop_autobsseog( EEG, [98], [98], 'bsscca', {'eigratio', [1000000]}, 'eog_fd', {'range',[1  4]});
    EEG = pop_autobsseog( EEG, [10], [10], 'sobi', {'eigratio', [1000000]}, 'eog_fd', {'range',[0  5]});
    alleeg(5) = EEG;
    
    % #### 3: Epoch extraction
    %epochs = extractEpochInterval(device);
    epochs = [0 5];
    EEG = pop_epoch( EEG, {'S769' 'S770'}, epochs, 'newname', ...
        'Alistair epochs', 'epochinfo', 'yes');
    alleeg(6) = EEG;
end


