function findCorrelationThreshold()
% Find thersholds for correlation coefficient
% After what correlation coefficient user start to feel that all channels are way too correlated
%The SpikeCleaner takes that threshold and sees if >80% of channels are as correlated as that.

% Here we run all the threshold values in a loop and match it with the user's curation, 
% to see at which threshold user seem to agree the most with the algorithm.

correlation=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.85,0.9,0.95,0.96,0.97,0.98,0.99];
correlationtsv=fullfile(pwd,'SpikeCleaner','cluster_Spikereasons.tsv');
baseThresholds = { ...
    'mode', 'lenient', ...
    'maxHW', 0.8, ...
    'minAmp', 50, ...
    'maxAmp', 2000, ...
    'minSlope', 150, ...
    'firingThreshold', 0.05, ...
    'acgallthreshold', 1.1, ...
    'acgallLabel', 'Noise', ...
    'corrThreshold', NaN ...
};

Correlation=fullfile(pwd,'Correlation');
if ~exist(Correlation,'dir')
    mkdir(Correlation);
end    
corrIdx = find(strcmp(baseThresholds,'corrThreshold')) + 1;
for i = 1:numel(correlation)
      
    baseThresholds{corrIdx} = correlation(i);     % set minHW for this run
    dz_classifyAllUnits(baseThresholds{:});

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
set(fig1, 'Position', [100 100 1100 700]);
set(fig1,'Units','pixels');
plot(thresholds_sorted, percentMatches_sorted, '-o'); 
xlabel('Max correlation threshold');
ylabel('(High correlation reason AND user=noise) in Percentage');
title('HighCorr reason vs user-noise matches across thresholds');
xLimits=xlim;
yLimits=ylim;
xPos = xLimits(1) + 0.05 * diff(xLimits);
yPos = yLimits(1) + 0.7 * diff(yLimits);

text(xPos, yPos, ...
    {'Choose the first value on x-axis, where y-axis value elbows and then plateaus.', ...
     'We are trying to find a coefficient of correlation at and after which user thinks all waveform channels are too correlated.', ...
     'You should see a rising curve that plateaus.'}, ...
      'FontSize', 10, ...
    'BackgroundColor', 'w', ...
    'EdgeColor', 'k', ...
    'Margin', 10);


set(fig1,'PaperPositionMode','auto');
print(fig1, fullfile(pwd,'Correlation', 'Correlation_NoiseThreshold'), '-dsvg');


%saving
saveas(fig1, fullfile(pwd,'Correlation', 'Correlation_NoiseThreshold.fig'));
saveas(fig1, fullfile(pwd,'Correlation', 'Correlation_NoiseThreshold.png'));

%% save as csv
corrTable = table( ...
    thresholds_sorted(:), ...
    percentMatches_sorted(:), ...
    'VariableNames', {'CorrelationThreshold','PercentNoiseMatches'});

writetable(corrTable, fullfile(pwd,'Correlation','Correlation_NoiseThreshold.csv'));
end




