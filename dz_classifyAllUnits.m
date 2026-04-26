function dz_classifyAllUnits(varargin)
% By Diksha Zutshi
% Function to run dz_Curate in the animal folder.

% Example:
% pipeline={'lowFiring','amplitude','halfWidth','slope','correlation','acgEmpty','acgAll','mua'};
% dz_classifyAllUnits('mode','strict','maxHW',0.8,'minAmp',50...,'pipeline'=pipeline)
 
    %look for all kilosort files in SpikeCleaner folder
    spikeCleaner=fullfile(pwd,'SpikeCleaner');
    if isfolder(spikeCleaner)
         % Detect .npy and .dat files
        npyFiles = dz_getAllExtFiles(spikeCleaner, 'npy', 0); % 1 = include subfolders
        datFiles = dz_getAllExtFiles(pwd, 'dat', 1);
    else
        error('Please run dz_runFirst.m to create SpikeCleaner Folder');
    end    
 
    [~,basename]=fileparts(pwd);
    datfilematch=[basename '.dat'];
   
    matchIdx=find(contains(datFiles,datfilematch),1,'first');
    if ~isempty(matchIdx)
        datfile=datFiles{matchIdx};
    else
        %
        error('Dat file isnt available: %s',datfilematch);      
    end    
    fname=datfile;
    
    % user can choose the number of channels to be considerer
    nLocalChannels = input('Enter number of local channels to consider: ');

    if isempty(nLocalChannels) || ~isscalar(nLocalChannels) || nLocalChannels < 1
        error('Please enter a valid positive integer for number of local channels.');
    end
    fprintf('Using local channels around max waveform channel: Example: 11 channels mean: Max channel + 10 closest channels to it.');
    fprintf('Enter %d how many channels would you like to consider according to your probe and shank geometry:.\n\n', nLocalChannels);
     % Parse name-value pair inputs
    p = inputParser;

    addParameter(p, 'mode', 'strict');              % ACG evaluation mode
    addParameter(p, 'maxHW', 0.8);                  % Half-width threshold in ms
    addParameter(p, 'minAmp', 50);                  % Minimum amplitude in uV
    addParameter(p, 'maxAmp', 500);                % Maximum amplitude in uV
    addParameter(p, 'minSlope', 150);               % Minimum slope in uV/ms
    addParameter(p, 'firingThreshold', 0.05);       % Minimum firing rate in Hz
    addParameter(p, 'acgallthreshold', 1);        % Threshold for all center bins vs shoulder
    addParameter(p, 'acgalllabel', 'Noise');        % Label for acgany threshold
    addParameter(p, 'corrThreshold', 0.95);         % Max correlation threshold
    addParameter(p, 'pipeline', ...
    {'lowFiring','correlation','amplitude','halfWidth','slope','acgEmpty','acgAll','mua'});



    parse(p, varargin{:});
    R = p.Results;

    % Convert parsed values into the cell format expected by dz_Curate
    thresholds = { ...
        R.mode, ...
        R.maxHW, ...
        R.minAmp, ...
        R.maxAmp, ...
        R.minSlope, ...
        R.firingThreshold, ...
        R.acgallthreshold, ...
        R.acgalllabel, ...
        R.corrThreshold ...     
    };        
    %Extracting pipeline after parsing
    pipeline = R.pipeline;
    

    %% Process .dat file

    fprintf('Processing file: %s\n', datfile);

    % Find the corresponding spike_clusters.npy file
    clusterfile = npyFiles(contains(npyFiles, 'spike_clusters'));
    clufile = clusterfile{1};  % Select the first match

    % Run curation
    dz_Curate(basename,fname, clufile,thresholds,pipeline,nLocalChannels)

    fprintf('Finished processing file: %s\n', fname);

end
