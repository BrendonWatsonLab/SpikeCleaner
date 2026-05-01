function dz_classifyAllUnits(varargin)
% By Diksha Zutshi
% Function to run dz_Curate in the animal folder.
% Run inside animal folder. Where you can see SpikeCleaner as a sub-folder

    %% Look for all Kilosort files in SpikeCleaner folder
    spikeCleaner = fullfile(pwd,'SpikeCleaner');

    if isfolder(spikeCleaner)
        npyFiles = dz_getAllExtFiles(spikeCleaner, 'npy', 0);
        datFiles = dz_getAllExtFiles(pwd, 'dat', 1);
    else
        error('Please run dz_runFirst.m to create SpikeCleaner folder');
    end

    [~, basename] = fileparts(pwd);
    datfilematch = [basename '.dat'];

    matchIdx = find(contains(datFiles, datfilematch), 1, 'first');

    if ~isempty(matchIdx)
        datfile = datFiles{matchIdx};
    else
        error('Dat file is not available: %s', datfilematch);
    end

    fname = datfile;

    %% Settings file path
    settingsFile = fullfile(spikeCleaner, sprintf('%s_SpikeCleaner_settings.json', basename));

    useSavedSettings = false;

    if exist(settingsFile, 'file')
        ansLoad = input('Existing SpikeCleaner settings found. Load them? y/n: ', 's');
        useSavedSettings = strcmpi(ansLoad, 'y');
    end

    %% Parser defaults
    p = inputParser;

    addParameter(p, 'mode', 'strict');
    addParameter(p, 'maxHW', 0.8);
    addParameter(p, 'minAmp', 50);
    addParameter(p, 'maxAmp', 500);
    addParameter(p, 'minSlope', 150);
    addParameter(p, 'firingThreshold', 0.05);
    addParameter(p, 'acgallthreshold', 1);
    addParameter(p, 'acgalllabel', 'Noise');
    addParameter(p, 'corrThreshold', 0.95);
    addParameter(p, 'pipeline', ...
        {'lowFiring','correlation','amplitude','halfWidth','slope','acgEmpty','acgAll','mua'});

    %% Load saved settings OR ask user
    if useSavedSettings
        oldSettings = jsondecode(fileread(settingsFile));

        varargin = { ...
            'mode', oldSettings.mode, ...
            'maxHW', oldSettings.maxHW, ...
            'minAmp', oldSettings.minAmp, ...
            'maxAmp', oldSettings.maxAmp, ...
            'minSlope', oldSettings.minSlope, ...
            'firingThreshold', oldSettings.firingThreshold, ...
            'acgallthreshold', oldSettings.acgallthreshold, ...
            'acgalllabel', oldSettings.acgalllabel, ...
            'corrThreshold', oldSettings.corrThreshold, ...
            'pipeline', cellstr(oldSettings.pipeline) ...
        };

    end

   

    %% Parse inputs
    parse(p, varargin{:});
    R = p.Results;

    %% Convert parsed values into threshold cell
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

    pipeline = R.pipeline;

    %% Save thresholds and pipeline for this run
    settings = struct();
    settings.mode = R.mode;
    settings.maxHW = R.maxHW;
    settings.minAmp = R.minAmp;
    settings.maxAmp = R.maxAmp;
    settings.minSlope = R.minSlope;
    settings.firingThreshold = R.firingThreshold;
    settings.acgallthreshold = R.acgallthreshold;
    settings.acgalllabel = R.acgalllabel;
    settings.corrThreshold = R.corrThreshold;
    settings.pipeline = R.pipeline;
    settings.dateRun = datestr(now);

    fid = fopen(settingsFile, 'w');

    if fid == -1
        error('Could not open settings file for writing: %s', settingsFile);
    end

    fprintf(fid, '%s', jsonencode(settings));
    fclose(fid);

    fprintf('Saved SpikeCleaner settings to:\n%s\n', settingsFile);

    %% Process .dat file
    fprintf('Processing file: %s\n', datfile);

    clusterfile = npyFiles(contains(npyFiles, 'spike_clusters'));

    if isempty(clusterfile)
        error('Could not find spike_clusters.npy in SpikeCleaner folder.');
    end

    clufile = clusterfile{1};

    dz_Curate(basename, fname, clufile, thresholds, pipeline);

    fprintf('Finished processing file: %s\n', fname);

end