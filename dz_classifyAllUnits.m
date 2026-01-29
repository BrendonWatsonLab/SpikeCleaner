function dz_classifyAllUnits(thresholds)
% By Diksha Zutshi
% Function to run dz_Curate on all .dat files in the current folder and subfolders.
% - mainDir: path to the code folder
% - thresholds: optional struct or cell array of thresholds
    
 
    %look for all kilosort files in SpikeCleaner folder
    spikeCleaner=fullfile(pwd,'SpikeCleaner');
    if isfolder(spikeCleaner)
         % Detect .npy and .dat files
        npyFiles = getAllExtFiles(spikeCleaner, 'npy', 0); % 1 = include subfolders
        datFiles = getAllExtFiles(pwd, 'dat', 1);
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
    
    %Default thresholds
    defaultThresholds = {...
        'lenient', ...      % ACG evaluation mode ('strict' or 'lenient')
        0.85, ...           % minHW: Half-width threshold in ms
        50, ...             % minAmp: Minimum amplitude in uV
        2000, ...           % maxAmp: Maximum amplitude in uV
        5e5, ...            % minSlope: Minimum slope (uV/s)
        0.05, ...           % firingThreshold: Minimum firing rate (Hz)
        0.8, ...            % acgallthreshold: Threshold for all center bins vs shoulder
        1.1 ...             % acgmaxthreshold: Threshold for any center bin vs shoulder
        0.95                %correlation 
   };

    % Use defaults if no thresholds provided
    if nargin < 2 || isempty(thresholds)
        thresholds = defaultThresholds;
    end

    %% Process .dat file

    fprintf('Processing file: %s\n', datfile);

    % Find the corresponding spike_clusters.npy file
    clusterfile = npyFiles(contains(npyFiles, 'spike_clusters'));
    clufile = clusterfile{1};  % Select the first match

    % Run curation
    dz_Curate(basename,fname, clufile,thresholds)

    fprintf('Finished processing file: %s\n', fname);

end
