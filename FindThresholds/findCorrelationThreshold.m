function findCorrelationThreshold()
% Find thersholds for correlation coefficient
% After what correlation coefficient user start to feel that all channels are way too correlated
%The SpikeCleaner takes that threshold and sees if >80% of channels are as correlated as that.

% Here we run all the threshold values in a loop and match it with the user's curation, 
% to see at which threshold user seem to agree the most with the algorithm.
correlationtsv=fullfile(pwd,'SpikeCleaner','cluster_Spikereasons.tsv');
correlation=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.85,0.9,0.95,0.96,0.97,0.98,0.99];
baseThresholds = { ...
    'lenient', ...   % ACG evaluation mode ('strict' or 'lenient')
    0.4, ...         % maxHW: Half-width threshold in ms (will be overwritten)
    50, ...          % minAmp: Minimum amplitude in uV
    2000, ...        % maxAmp: Maximum amplitude in uV
    150, ...         % minSlope: Minimum slope (uV/s)
    0.05, ...        % firingThreshold: Minimum firing rate (Hz)
    0.9, ...         % acgallthreshold: Threshold for all center bins vs shoulder
    Noise, ...         % label for acgall: Noise or MUA:::from results from findACGallThreshold()
    1.1 , ...         % acgmaxthreshold: Threshold for any center bin vs shoulder
    Noise, ...         % label for acgany: Noise or MUA::: from results from findACGanyThreshold()..          
    NaN              % correlation
};




Correlation=fullfile(pwd,'Correlation');
if ~exist(Correlation,'dir')
    mkdir(Correlation);
end    
for i = 1:numel(correlation)
      
    baseThresholds{9} = correlation(i);     % set minHW for this run
    dz_classifyAllUnits(baseThresholds);

    %renaming the tsv spikecleaner outputs to halfwidth_spikecleaner
    dstFile = fullfile(Correlation, sprintf('Correlation_%0.2f.tsv', correlation(i)));
    movefile(correlationtsv, dstFile);
end
%% Plotting

tsvfiles=dir(fullfile(pwd,'Correlation','*.tsv'));
files = dir(fullfile(pwd,'SpikeCleaner', 'cluster_group*.tsv'));
filePath = fullfile(files(1).folder, files(1).name);
T1=readtable(filePath, ...
        'FileType','text', ...
        'Delimiter','\t');
thresholds = nan(length(tsvfiles),1);  
for i=1:length(tsvfiles)
    algo = fullfile(tsvfiles(i).folder, tsvfiles(i).name);
    tok = regexp(tsvfiles(i).name, '^Correlation_(\d*\.?\d+)\.tsv$', 'tokens', 'once');

    thresholds(i) = str2double(tok{1});
    
    T = readtable(algo, ...
        'FileType','text', ...
        'Delimiter','\t');
    
    % Cleaning up strings
    T.Reasons = strtrim(lower(T.Reasons));
    % Merge tables on cluster_id
    merged = innerjoin(T, T1, 'Keys', 'cluster_id');      
     % ---- High correlation reason rows (in the *merged* table)
    idxHighCorr = startsWith(merged.Reasons, "high correlation");  % already lowercased
    
    vars = merged.Properties.VariableNames;
    groupCol = vars(contains(lower(vars),'group'));
    userLabels = merged.(groupCol{1});
    % ---- user thinks noise (in merged.group)
    idxUserNoise = (lower(strtrim(string(userLabels))) == "noise");

    % ---- match: high corr reason & user noise
    highCorrNoiseMatches(i) = sum(idxHighCorr == idxUserNoise);   
    
    % total units
    totalUnits = height(merged);

    % percentage
    percentMatches(i) = 100 * highCorrNoiseMatches(i) / totalUnits;
end
% 
% remove skipped
valid = ~isnan(thresholds);
thresholds = thresholds(valid);
percentMatches = percentMatches(valid);
[thresholds_sorted, idx] = sort(thresholds, 'ascend');
percentMatches_sorted = percentMatches(idx);


% % plot
fig1=figure;
plot(thresholds_sorted, percentMatches_sorted, '-o'); 
xlabel('Max correlation threshold');
ylabel('(High correlation reason AND user=noise) in Percentage');
title('HighCorr reason vs user-noise matches across thresholds');
xLimits=xlim;
yLimits=ylim;
text(xLimits(2)*0.3,yLimits(2)*0.3, {'Choose the first value on x-axis, where y-axis value elbows and then platues.', ...
               'We are trying to find a  coeffient of correlation at and after which user thinks all waveform channels are too correlated.', ...
               ' You should see a rising curve that platues.'}, ...
     'FontSize', 10);
%saving
saveas(fig1, fullfile(pwd,'Correlation', 'Correlation_NoiseThreshold.png'));


% % 
end




