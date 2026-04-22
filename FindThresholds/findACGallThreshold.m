function findACGallThreshold()
% Find thersholds for ACGall: above that percentage of fill in ACG user thinks ACG is too full to be neuronal.

% We expect to see a rising curve which plateaus: matches increasing as the threshold rising, 
%as it's a one sided threshold

% Here we run all the threshold values in a loop and match it with the user's curation, 
% to see at which threshold user seem to agree the most with the algorithm.

% acgallsv=fullfile(pwd,'SpikeCleaner','cluster_Spikereasons.tsv');
% % % %% run classifyAllUnits with different ACG thresholds
% acgmax1=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5];
% 
% baseThresholds = { ...
%     'mode', 'lenient', ...
%     'maxHW', 0.8, ...
%     'minAmp', 50, ...
%     'maxAmp', 500, ...
%     'minSlope', 150, ...
%     'firingThreshold', 0.05, ...
%     'acgallthreshold', NaN, ...
%     'acgalllabel', 'Noise', ...
%     'corrThreshold', 0.95 ...
% };
% 
%   
% 
% ACGall=fullfile(pwd,'ACGall');
% if ~exist(ACGall,'dir')
%     mkdir(ACGall);
% end    
% AcgallIdx = find(strcmp(baseThresholds,'acgallthreshold')) + 1;
% 
% for i = 1:numel(acgmax1)
%     baseThresholds{AcgallIdx} = acgmax1(i);     
%     dz_classifyAllUnits(baseThresholds{:});
%     
%     dstFile = fullfile(ACGall, sprintf('ACGall_%0.2f.tsv', acgmax1(i)));
%     movefile(acgallsv, dstFile);
% end


%% Plotting
tsvfiles=dir(fullfile(pwd,'ACGall','*.tsv'));
files = dir(fullfile(pwd,'SpikeCleaner', 'cluster_group*.tsv'));
filePath = fullfile(files(1).folder, files(1).name);
T1=readtable(filePath, ...
        'FileType','text', ...
        'Delimiter','\t');
thresholds = nan(length(tsvfiles),1);   
 
for i=1:length(tsvfiles)
    algo = fullfile(tsvfiles(i).folder, tsvfiles(i).name);
    tok = regexp(tsvfiles(i).name, '^ACGall_(\d*\.?\d+)\.tsv$', 'tokens', 'once');

    thresholds(i) = str2double(tok{1});
    T = readtable(algo, ...
        'FileType','text', ...
        'Delimiter','\t');
    
      % Cleaning up strings
    T.Reasons = strtrim(lower(T.Reasons));
    % Merge tables on cluster_id
    merged = innerjoin(T, T1, 'Keys', 'cluster_id');
        
     % ---- High correlation reason rows (in the *merged* table)
    idxHighCorr = startsWith(merged.Reasons, "centerbins proportion greater");  
    
     % renaming the merged user label
    vars = merged.Properties.VariableNames;
    groupCol = vars(contains(lower(vars),'group'));
    userLabels = merged.(groupCol{1});
    %% user thinks such cases are noise
    % ---- user thinks noise (in merged.group)
    idxUserNoise = (lower(strtrim(string(userLabels))) == "noise");    
    % ---- match: high corr reason & user noise
    highCorrNoiseMatches(i) = sum(idxHighCorr == idxUserNoise);       
    % total units
    totalUnits = height(merged);
    % percentage
    percentMatchesnoise(i) = 100 * highCorrNoiseMatches(i) / totalUnits;        
    
    %% user thinks such caes are mua
    % ---- user thinks noise (in merged.group)
    idxUserNoise1 = (lower(strtrim(string(userLabels))) == "mua");   
    % ---- match: high corr reason & user noise
    highCorrNoiseMatches1(i) = sum(idxHighCorr == idxUserNoise1);   
    % percentage
    percentMatchesmua(i) = 100 * highCorrNoiseMatches1(i) / totalUnits;    
        
end

betterLabel = strings(length(thresholds),1);

for i = 1:length(thresholds)
    if percentMatchesnoise(i) > percentMatchesmua(i)
        betterLabel(i) = "noise";
    elseif percentMatchesmua(i) > percentMatchesnoise(i)
        betterLabel(i) = "mua";
    else
        betterLabel(i) = "tie";
    end
end

resultTable = table(thresholds(:), percentMatchesnoise(:), percentMatchesmua(:), betterLabel(:), ...
    'VariableNames', {'Threshold','NoiseMatchPct','MuaMatchPct','BetterLabel'});

disp(resultTable)

%% choose the rising curve between mua and noise labels

%noise
% remove skipped
valid = ~isnan(thresholds);
thresholds = thresholds(valid);
percentMatches = percentMatchesnoise(valid);
[thresholds_sorted, idx] = sort(thresholds, 'ascend');
percentMatches_sorted1 = percentMatches(idx);
% 
% % % plot
% figure;
% plot(thresholds_sorted, percentMatches_sorted1, '-o'); grid on;
% xlabel('ACG all bins threshold');
% ylabel('Total matches for Noise labels(Only counting all ACG Bins being higher than the thershold)');
% title('Max allowed percentage in all centerbins required for users to think its still neuronal: Otherwise Noise');
% grid on;

%mua
% remove skipped
valid = ~isnan(thresholds);
thresholds = thresholds(valid);
percentMatches = percentMatchesmua(valid);
[thresholds_sorted, idx] = sort(thresholds, 'ascend');
percentMatches_sorted2 = percentMatches(idx);
% 
% % % plot
% figure;
% plot(thresholds_sorted, percentMatches_sorted2, '-o'); grid on;
% xlabel('ACG all bins threshold');
% ylabel('Total matches for MUA labels(Only counting all ACG Bins being higher than the thershold)');
% title('Max allowed percentage in all centerbins required for users to think its still a Single Unit-: Otherwise MUA');
% grid on;


% both together
fig1 = figure;
set(fig1, 'Position', [100 100 1100 700]);
set(fig1,'Units','pixels');
plot(thresholds_sorted, percentMatches_sorted1, '-o'); hold on; grid on;
plot(thresholds_sorted, percentMatches_sorted2, '-s'); 

xlabel('ACG all bins threshold');
ylabel('Total matches for Noise/MUA labels');
title({'Max allowed percentage in all center bins required for users to think it is still neuronal', ...
       'Otherwise classified as Noise / MUA'});

legend('Noise', 'MUA', 'Location', 'best');

xLimits = xlim;
yLimits = ylim;

xPos = xLimits(1) + 0.05 * diff(xLimits);
yPos = yLimits(1) + 0.70 * diff(yLimits);

text(xPos, yPos, ...
    {'Choose the label of the rising curve for this criterion and plug it into defaultThresholds in dzclassifyAllUnits.', ...
     'Choose the first x-value where the y-value elbows and then plateaus.'}, ...
    'FontSize', 10, ...
    'BackgroundColor', 'w', ...
    'EdgeColor', 'k', ...
    'Margin', 10);

set(gca,'LooseInset',max(get(gca,'TightInset'),0.02));
set(fig1,'PaperPositionMode','auto');


print(fig1, fullfile(pwd,'ACGall','ACGall_Noise_vs_MUA'), '-dsvg');









%saving
% saveas(fig1, fullfile(pwd,'ACGall','ACGall_Noise_vs_MUA.svg'));

% saveas(fig1, fullfile(pwd,'ACGall', 'ACGall_Noise_vs_MUA.fig'));

%% save as csv
acgallTable = table( ...
    thresholds_sorted(:), ...
    percentMatches_sorted1(:), ...
    percentMatches_sorted2(:), ...
    'VariableNames', {'ACGallThreshold','PercentNoiseMatches','PercentMUAMatches'});

writetable(acgallTable, fullfile(pwd,'ACGall','ACGall_NoiseThreshold.csv'));




end