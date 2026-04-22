function findACGanyThreshold()
% Find thersholds for ACGany: above that percentage of fill in ACG user thinks ACG is too full to be neuronal

% Here we run all the threshold values in a loop and match it with the user's curation, 
% to see at which threshold user seem to agree the most with the algorithm.

acganysv=fullfile(pwd,'SpikeCleaner','cluster_Spikereasons.tsv');
% % %% run classifyAllUnits with different ACG thresholds
acgmax2=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5];

baseThresholds = { ...
    'mode', 'lenient', ...
    'maxHW', 0.7, ...
    'minAmp', 50, ...
    'maxAmp', 2000, ...
    'minSlope', 150, ...
    'firingThreshold', 0.05, ...
    'acganythreshold', NaN, ...
    'acganyLabel', 'Noise', ...
    'corrThreshold', 0.95 ...
};

ACGany=fullfile(pwd,'ACGany');
if ~exist(ACGany,'dir')
    mkdir(ACGany);
end   
AcganyIdx = find(strcmp(baseThresholds,'acganythreshold')) + 1;

for i = 1:numel(acgmax2)
    
    baseThresholds{AcganyIdx} = acgmax2(i);     % set minHW for this run
    dz_classifyAllUnits(baseThresholds{:});
    
    dstFile = fullfile(ACGany, sprintf('ACGany_%0.2f.tsv', acgmax2(i)));
    movefile(acganysv, dstFile);
end


%% Plotting
tsvfiles=dir(fullfile(pwd,'ACGany','*.tsv'));
files = dir(fullfile(pwd,'SpikeCleaner', 'cluster_group*.tsv'));
filePath = fullfile(files(1).folder, files(1).name);
T1=readtable(filePath, ...
        'FileType','text', ...
        'Delimiter','\t');
thresholds = nan(length(tsvfiles),1);   
for i=1:length(tsvfiles)
    algo = fullfile(tsvfiles(i).folder, tsvfiles(i).name);
    tok = regexp(tsvfiles(i).name, '^ACGany_(\d*\.?\d+)\.tsv$', 'tokens', 'once');

    thresholds(i) = str2double(tok{1});
    T = readtable(algo, ...
        'FileType','text', ...
        'Delimiter','\t');
    
      % Cleaning up strings
    T.Reasons = strtrim(lower(T.Reasons));
    % Merge tables on cluster_id
    merged = innerjoin(T, T1, 'Keys', 'cluster_id');
        
     % ---- High correlation reason rows (in the *merged* table)
    idxHighCorr = startsWith(merged.Reasons, "atleast one");  % already lowercased

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
disp('Choose the label for the chosen Threshold value for this criteria');

%% choose the rising curve between mua and noise labels
%noise
% remove skipped
valid = ~isnan(thresholds);
thresholds = thresholds(valid);
percentMatches = percentMatchesnoise(valid);
[thresholds_sorted, idx] = sort(thresholds, 'ascend');
percentMatches_sorted1 = percentMatches(idx);


plot(thresholds_sorted, percentMatches_sorted1, '-o'); grid on;
xlabel('ACG any bins threshold');
ylabel('Total matches for Noise labels(Only counting any ACG Bins being higher than the thershold)');
title('Max allowed percentage in any centerbins required for users to think its still neuronal: Otherwise Noise');
grid on;


%noise
% remove skipped
valid = ~isnan(thresholds);
thresholds = thresholds(valid);
percentMatches = percentMatchesmua(valid);
[thresholds_sorted, idx] = sort(thresholds, 'ascend');
percentMatches_sorted2 = percentMatches(idx);


plot(thresholds_sorted, percentMatches_sorted2, '-o'); grid on;
xlabel('ACG any bins threshold');
ylabel('Total matches for MUA labels(Only counting any ACG Bins being higher than the thershold)');
title('Max allowed percentage in any centerbins required for users to think its still a Single Unit-: Otherwise MUA');
grid on;


% both together
fig1 = figure;
plot(thresholds_sorted, percentMatches_sorted1, '-o'); grid on;
hold on;
plot(thresholds_sorted, percentMatches_sorted2, '-s'); grid on;
xlabel('ACG any bins threshold');
ylabel('Total matches for MUA/Noise labels(Only counting any ACG Bins being higher than the thershold)');
title('Max allowed percentage in any centerbins required for users to think its still a Single Unit-: Otherwise MUA/Noise');
legend('Noise', 'MUA', 'Location', 'best');

text(0.6, 80, {'Choose the label of the rising curve for this criteria and plug it in defaultThresholds in dz_ classifyAllUnits .', ...
               'Choose the first value on x-axis, where y-axis value elbows and then platues.'}, ...
     'FontSize', 10);
%grid on;

%saving
saveas(fig1, fullfile(pwd,'ACGany', 'ACGany_Noise_vs_MUA.png'));





end