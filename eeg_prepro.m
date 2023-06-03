% EEG Preprocessing Pipeline Step 1
% Matt Kmiecik
% Started 03 June 2023

workspace_prep % Prepares workspace (see /src)

% Initializes subjects for batch processing (if applicable)
ss = string({RAW{2:size(RAW,1),1}});

i=1; % for testing purposes

% Preprocessing ----
for i = 1:length(ss)

    % Creating variables ----
    this_ss = NUM(i);
    this_ss_path = dir(fullfile(data_dir, strcat(num2str(this_ss), '.bdf')));
    this_ss_name = this_ss_path.name;

    % Loads in raw data using biosemi ----
    EEG = pop_biosig(...
        fullfile(this_ss_path.folder, this_ss_name),...
        'ref',[1] ,...
        'refoptions',{'keepref','on'}...
        );
    
    
    % Remove externals that are not being used ----
    EEG = pop_select( EEG, 'rmchannel',{'EXG3','EXG4','EXG5','EXG6','EXG7','EXG8'});
    
    % Configuring channel locations ----
    % loads in ELP
    eloc = readlocs( this_elp ); % reads in elp chan locations
    EEG.chanlocs = eloc(4:69); % adds these except for fiducials
    % sets A1 as ref because it was chosen upon import
    EEG = pop_chanedit(EEG, 'setref', {'1:66' 'A1'}); 
    
    % Downsamples to 256Hz ----
    EEG = pop_resample(EEG, 256);
    
    % Re-references data to linked mastoids
    EEG = pop_reref( EEG, [24 61] ,'keepref','on');
    
    % Removing DC offset by subtracting the mean signal from each electrode
    EEG = pop_rmbase(EEG, [], []);
    
    % Highpass filter at 1Hz 
    % -6dB @ 1Hz, 425 point highpass, 2Hz transition band width
    EEG = pop_eegfiltnew(EEG, 'locutoff', 2, 'plotfreqz', 0);

    % Cleanline ----
    % Removing electrical line noise @ 60 Hz
    EEG = pop_cleanline(EEG, 'bandwidth', 2, 'chanlist', [1:EEG.nbchan],...
        'computepower', 1, 'linefreqs', 60, 'normSpectrum', 0, ...
        'p', 0.01, 'pad', 2, 'plotfigures' , 0, 'scanforlines', 1, ...
        'sigtype', 'Channels', 'tau', 100, 'verb', 1, 'winsize', ...
        4,'winstep',1);
    
    % Renames dataset ----
    dataset_name = strcat(num2str(this_ss), '-prepro');
    EEG = pop_editset(EEG, 'setname', dataset_name, 'run', []);

    % Saves out preprocessed data for inspection ----
    EEG = pop_saveset(EEG, 'filename', dataset_name, 'filepath', output_dir);
    
    eeglab redraw % redraws to GUI for convenience

end