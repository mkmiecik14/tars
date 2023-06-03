% EEG ICA Step 2 (after visual inspection)
% Matt Kmiecik
% Started 20 July 2022

workspace_prep % Prepares workspace

% Preprocessing ----
num_iters = size(NUM, 1); % number of loop iterations/participants
i=1; % for testing purposes

for i = 1:num_iters
    
    % Creating variables ----
    this_ss = NUM(i);
    this_ss_path = dir(fullfile(output_dir, strcat(num2str(this_ss), '*-prepro.set')));
    this_ss_name = this_ss_path.name;
    
    % Checks to see if the a visually inspected and rejected file exists
    % (i.e., if sections of the EEG were rejected due to noise)
    this_ss_vis_rej = strcat(num2str(this_ss), '-prepro-vis-rej.set');
    if isfile(fullfile(output_dir, this_ss_vis_rej))
     % File exists...load in the visually inspected and rejected file
     % Loads in raw data using EEGLAB ----
    EEG = pop_loadset('filename',this_ss_vis_rej,'filepath', this_ss_path.folder);
    else
     % File does not exist...load standard file
    % Loads in raw data using EEGLAB ----
    EEG = pop_loadset('filename',this_ss_name,'filepath', this_ss_path.folder);
    end
    
    % If there were any bad channels identified during inspection, they are
    % imported here:
    interpchans = str2num(RAW{i+1,2}); % gathers bad channels from excel
    
    % Interpolates bad channels if bad channels were identified
    if sum(size(interpchans, 2)) > 0 % checks to see if there are bad channels
            disp('Bad channel(s) detected...interpolating channels...');
            EEG = pop_interp(EEG, interpchans, 'spherical');
    else
            disp('No bad channels detected...')
    end
    
    % Extracts resting state blocks ----
    EEG = pop_rmdat( EEG, {'102'},[0 60] ,0);
    
    
    % ICA decomposition ----
    % computes rank from EEG.nbchan( == 30) - n (interpolated chans)
    this_rank = EEG.nbchan - 1 - length(interpchans); % minus 1 due to ref
    EEG = pop_runica(...
        EEG,...
        'icatype', 'runica',... 
        'extended',1,...
        'interrupt','on',...
        'pca',this_rank...
    );
    
    % renames dataset
    dataset_name = strcat('rs-', visit_name, '-', num2str(this_ss), '-ica');
    EEG = pop_editset(EEG, 'setname', dataset_name, 'run', []);
    
    % Saves out data
    outname = strcat(dataset_name, '.set'); % save out subject name
    EEG = pop_saveset(EEG, 'filename', outname, 'filepath', outpath);
    
end

eeglab redraw