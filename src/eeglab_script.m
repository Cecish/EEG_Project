% Instead of using the EEGLAB GUI
% Params:
%   - path: path of the raw EEG data file
%   - name: name of the raw EEG data file
%   - name_dataset: name of the initial data set to create
%   - highpass: highpass filter
%   - notch: notch filter
%   - epochs: window interval for epoch extraction
% Return: a alleeg structure similar to the ALLEEG variable produced by
% eeglab toolbox when using the GUI
function alleeg = eeglab_func(path, name, name_dataset, highpass, notch, epochs)
    eegh;
    
    % #### 0: Load the data file into EEGLAB toolbox
    EEG = pop_loadbv(path, name);
    EEG.setname = name_dataset;
    
    % #### 1: Edit channels location
    EEG = pop_chanedit(EEG, 'lookup');
    alleeg(1) = EEG;
    
    % #### 2: Preprocessing
    % 2.1 Re-reference (average) ---------> 23.8934%
    EEG = pop_reref( EEG, []);
    EEG.setname = 'Alistair re-referencing average';
    alleeg(2) = EEG;
    
    % 2.2 High pass filter (1Hz) ---------> 20.9124%
    EEG = pop_eegfiltnew(EEG, [], highpass, 424, true, [], 0);
    EEG.setname = 'Alistair highpass filter';
    alleeg(3) = EEG;
    
    % 2.3 Notch filter ------------> 20.8672%
    EEG = pop_eegfiltnew(EEG, notch(1), notch(2), 212, 1, [], 0);
    EEG.setname = 'Alistair notch filter';
    alleeg(4) = EEG;
    
    EEG = pop_autobsseog( EEG, [98], [98], 'bsscca', {'eigratio', [1000000]}, 'eog_fd', {'range',[1  4]});
    alleeg(5) = EEG;
    
    % #### X: Epoch extraction
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

