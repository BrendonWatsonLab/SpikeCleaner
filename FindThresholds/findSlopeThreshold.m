function findSlopeThreshold()
% Find thersholds for slope: under that point user thinks waveform isnt
% neuronal because the slope is too low, waveform is too slow

% Here we run all the threshold values in a loop and match it with the user's curation, 
% to see at which threshold user seem to agree the most with the algorithm.
% 
slopesstsv=fullfile(pwd,'SpikeCleaner','cluster_Spikereasons.tsv');
slopes=[0.1,0.5,1,5,10,50,100,150,200,250,300,350,400 500,600,700,800,900,1000,1500,2000];
baseThresholds = { ...
    'mode', 'lenient', ...
    'maxHW', 0.8, ...
    'minAmp', 50, ...
    'maxAmp', 1000, ...
    'minSlope', NaN, ...
    'firingThreshold', 0.05, ...
    'acgallthreshold', 1.1, ...
    'acgallLabel', 'Noise', ...
    'corrThreshold', 0.95 ...
};
Slope=fullfile(pwd,'Slope');
if ~exist(Slope,'dir')
    mkdir(Slope);
end    
slopeIdx = find(strcmp(baseThresholds,'minSlope')) + 1;
for i = 1:numel(slopes)
     
    baseThresholds{slopeIdx} = slopes(i);     % set minHW for this run
    dz_classifyAllUnits(baseThresholds{:});

    %renaming the tsv spikecleaner outputs to halfwidth_spikecleaner
    dstFile = fullfile(Slope, sprintf('slope_%0.2f.tsv', slopes(i)));
    movefile(slopesstsv, dstFile);
end


%% Plotting
tsvfiles=dir(fullfile(pwd,'Slope','*.tsv'));
files = dir(fullfile(pwd,'SpikeCleaner', 'cluster_group*.tsv'));
filePath = fullfile(files(1).folder, files(1).name);
T1=readtable(filePath, ...
        'FileType','text', ...
        'Delimiter','\t');
thresholds = nan(length(tsvfiles),1);   

for i=1:length(tsvfiles)
    algo = fullfile(tsvfiles(i).folder, tsvfiles(i).name);
    tok = regexp(tsvfiles(i).name, '^slope_(\d*\.?\d+)\.tsv$', 'tokens', 'once');

    thresholds(i) = str2double(tok{1});
    T = readtable(algo, ...
        'FileType','text', ...
        'Delimiter','\t');
    
      % Cleaning up strings
    T.Reasons = strtrim(lower(T.Reasons));
    % Merge tables on cluster_id
    merged = innerjoin(T, T1, 'Keys', 'cluster_id');
        
     % ---- High correlation reason rows (in the *merged* table)
    idxHighCorr = startsWith(merged.Reasons, "slope too");  % already lowercased

     % renaming the merged user label
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

% remove skipped
valid = ~isnan(thresholds);
thresholds = thresholds(valid);
percentMatches = percentMatches(valid);
[thresholds_sorted, idx] = sort(thresholds, 'ascend');
percentMatches_sorted = percentMatches(idx);


% % plot
fig1 = figure;
set(fig1, 'Position', [100 100 1100 700]);
set(fig1,'Units','pixels');
plot(thresholds_sorted, percentMatches_sorted, '-o');
xlabel('Slope threshold-uV/ms');
ylabel('Total matches for Noise labels(Slope)-Percetage Matches');
title('Min Slope required for users to think a spike is neuronal');
xLimits=xlim;
yLimits=ylim;
xPos = xLimits(1) + 0.05 * diff(xLimits);
yPos = yLimits(1) + 0.70 * diff(yLimits);

text(xPos, yPos, ...
    {'Choose the last value on x-axis, where y-axis value elbows and then falls.', ...
               'We are trying to find a waveform slope(uV/ms) lower than which a spike appears non physiological, too slow to progress.', ...
               'You should see a falling curve because matches would decrease if we choose a bigger minimum slope ',...
               'thershold as those spikes would still look physiological.'}, ...
    'FontSize', 10, ...
    'BackgroundColor', 'w', ...
    'EdgeColor', 'k', ...
    'Margin', 10);

set(fig1,'PaperPositionMode','auto');
print(fig1, fullfile(pwd,'Slope','Slope_NoiseThreshold'), '-dsvg');

%saving
saveas(fig1, fullfile(pwd,'Slope', 'Slope_NoiseThreshold.png'));
saveas(fig1, fullfile(pwd,'Slope', 'Slope_NoiseThreshold.fig'));
%% save as csv
slopeTable = table( ...
    thresholds_sorted(:), ...
    percentMatches_sorted(:), ...
    'VariableNames', {'SlopeThreshold','PercentNoiseMatches'});

writetable(slopeTable, fullfile(pwd,'Slope','Slope_NoiseThreshold.csv'));


end