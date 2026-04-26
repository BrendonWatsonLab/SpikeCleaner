function findAmplitudeThreshold()
% 
amplitudes = [100,200,300,400,500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,11000,12000,13000,14000,15000];
amplitudetsv=fullfile(pwd,'SpikeCleaner','cluster_Spikereasons.tsv');
baseThresholds = { ...
    'mode', 'lenient', ...
    'maxHW', 0.8, ...
    'minAmp', 50, ...
    'maxAmp', NaN, ...
    'minSlope', 150, ...
    'firingThreshold', 0.05, ...
    'acgallthreshold', 1.1, ...
    'acgallLabel', 'Noise', ...
    'corrThreshold', 0.95 ...
};


Amplitude=fullfile(pwd,'Amplitude');
if ~exist(Amplitude,'dir')
    mkdir(Amplitude);
end   
ampIdx = find(strcmp(baseThresholds,'maxAmp')) + 1;
for i = 1:numel(amplitudes)
      
    baseThresholds{ampIdx} = amplitudes(i);     % set minHW for this run
    dz_classifyAllUnits(baseThresholds{:});

    %renaming the tsv spikecleaner outputs to halfwidth_spikecleaner
    dstFile = fullfile(Amplitude, sprintf('Amplitude_%0.2f.tsv', amplitudes(i)));
    movefile(amplitudetsv, dstFile);
end

%% Plotting

tsvfiles=dir(fullfile(pwd,'Amplitude','*.tsv'));
files = dir(fullfile(pwd,'SpikeCleaner', 'cluster_group*.tsv'));
filePath = fullfile(files(1).folder, files(1).name);
T1=readtable(filePath, ...
        'FileType','text', ...
        'Delimiter','\t');
thresholds = nan(length(tsvfiles),1);  
for i=1:length(tsvfiles)
    algo = fullfile(tsvfiles(i).folder, tsvfiles(i).name);
    tok = regexp(tsvfiles(i).name, '^Amplitude_(\d*\.?\d+)\.tsv$', 'tokens', 'once');

    thresholds(i) = str2double(tok{1});
    
    T = readtable(algo, ...
        'FileType','text', ...
        'Delimiter','\t');
    
    % Cleaning up strings
    T.Reasons = strtrim(lower(T.Reasons));
    % Merge tables on cluster_id
    merged = innerjoin(T, T1, 'Keys', 'cluster_id');      
     % ---- High correlation reason rows (in the *merged* table)
    idxHighCorr = startsWith(merged.Reasons, "amplitude is too low");  % already lowercased
    
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
xlabel('Max amplitude threshold');
ylabel('(High amplitude reason AND user=noise) in Percentage');
title('Too high of amplitude vs user-noise matches across thresholds');
xLimits=xlim;
yLimits=ylim;
xPos = xLimits(1) + 0.05 * diff(xLimits);
yPos = yLimits(1) + 0.70 * diff(yLimits);

text(xPos, yPos, ...
   {'Choose the last value on x-axis, where y-axis value elbows and then falls.', ...
               'We are trying to find a threshold for amplitude at and after user thinks the difference in absolute', ... 
               'amplitude between channels is too much for it to be neuronal', ...
               'You should see a rising curve because matches would incraese if we choose a bigger', ...
               'maximum amplitude difference thershold as we will be missing labelling many units that user think are noisy.'}, ...
    'FontSize', 10, ...
    'BackgroundColor', 'w', ...
    'EdgeColor', 'k', ...
    'Margin', 10);

set(fig1,'PaperPositionMode','auto');
print(fig1, fullfile(pwd,'Amplitude', 'Amplitude_NoiseThreshold'), '-dsvg');


 
 
 
 
 
 
%saving
saveas(fig1, fullfile(pwd,'Amplitude', 'Amplitude_NoiseThreshold.png'));
end