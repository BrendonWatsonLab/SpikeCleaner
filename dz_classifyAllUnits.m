function dz_classifyAllUnits(thresholds)
% By Diksha Zutshi
% Function to run dz_Curate on all .dat files in the current folder and subfolders.
% - mainDir: path to the code folder
% - thresholds: optional struct or cell array of thresholds
    
    %addpath(genpath("D:\Whatson Lab\Klustakwik\Matlab\SpikeCleaner\HelperFunctions"));
    % Detect .npy and .dat files
    npyFiles = getAllExtFiles(pwd, 'npy', 1); % 1 = include subfolders
    datFiles = getAllExtFiles(pwd, 'dat', 1);
    basename = cellfun(@bz_BasenameFromBasepath, datFiles, 'UniformOutput', false);

    %Default thresholds
    defaultThresholds = {...
        'strict', ...       % ACG evaluation mode ('strict' or 'lenient')
        0.85, ...           % minHW: Half-width threshold in ms
        50, ...             % minAmp: Minimum amplitude in uV
        2000, ...           % maxAmp: Maximum amplitude in uV
        5e5, ...            % minSlope: Minimum slope (uV/s)
        0.05, ...           % firingThreshold: Minimum firing rate (Hz)
        0.8, ...            % acgallthreshold: Threshold for all center bins vs shoulder
        1.1 ...             % acgmaxthreshold: Threshold for any center bin vs shoulder
    };

    % Use defaults if no thresholds provided
    if nargin < 2 || isempty(thresholds)
        thresholds = defaultThresholds;
    end

    % Process .dat file
    for i = 1:length(datFiles)
        fname = datFiles{i};

        if exist(fname, 'file')
            try
                fprintf('Processing file: %s\n', fname);

                % Find the corresponding spike_clusters.npy file
                clusterfile = npyFiles(contains(npyFiles, 'spike_clusters'));
                if isempty(clusterfile)
                    warning('No spike_clusters.npy file found for %s. Skipping.', fname);
                    continue;
                end
                clufile = clusterfile{1};  % Select the first match

                % Run curation
                dz_Curate(basename,fname, clufile,thresholds)

                fprintf('Finished processing file: %s\n', fname);
            catch ME
                fprintf('Error processing file: %s\n', fname);
                fprintf('Error: %s\n', ME.message);
            end
        end
    end
end
