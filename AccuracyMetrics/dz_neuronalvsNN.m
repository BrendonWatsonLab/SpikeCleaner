function dz_neuronalvsNN(curdir)
if nargin < 1 || isempty(curdir)
    curdir=pwd;
end    
%%%%%%%%%this code is for comparing good vs rest: SU Vs Rest
%%this will export a metrics.tsv for match or no match
username=strtrim(lower(input('Enter the user name for the PHY display:','s')));
group1 =readtable("cluster_SpikeCleaner.tsv", 'FileType', 'text', 'Delimiter', '\t');%spike cleaner outputs
userlabels='cluster_group.tsv';
renameduser=sprintf('cluster_group_%s.tsv',username);
group= readtable(renameduser, 'FileType', 'text', 'Delimiter', '\t');
 % Rename header "group" to "username_group"
oldVar = 'group';
newVar = sprintf('%s_group', username);
if ~isfile(renameduser)    
    if isfile(userlabels)
        movefile(userlabels,renameduser); %renames the original tsv file
        fprintf('renamed the user labels %s to %s \n',userlabels,renameduser);
        %rename header name        
        idx = strcmp(group.Properties.VariableNames, oldVar);
        if any(idx)
            group.Properties.VariableNames{idx} = newVar;
        else
            warning('Column "%s" not found.', oldVar);
        end

        % Write to new file
        writetable(group, renameduser, 'FileType', 'text', 'Delimiter', '\t');

        fprintf('Created %s with updated column header.\n', renameduser);
    else
        warning('Cant find cluster_group.tsv, please manually curate using PHY');
    end    
end
 
% Cleaning up strings
group1.SpikeCleaner = strtrim(lower(group1.SpikeCleaner));

% Merge tables on cluster_id
merged = innerjoin(group, group1, 'Keys', 'cluster_id');

 %: "neuronal" if good or mua
true_label = ismember(merged.(newVar), ["good","mua"]);

% pred_label: 1 if not noise/lowrate, 0 otherwise
non_neuronal_algo = contains(merged.SpikeCleaner, "noise") | ...
                    contains(merged.SpikeCleaner, "lowrate");
pred_label = ~non_neuronal_algo;

% Initialize fp_fn column
merged.fp_fn = repmat("NotMatch", height(merged), 1);
mask = (true_label == pred_label);   % same decision
merged.fp_fn(mask) = "Match";


% Confusion matrix values
TP = sum(true_label == 1 & pred_label == 1);  % Correctly identified neuronal
FP = sum(true_label == 0 & pred_label == 1);  % Non-neuronal incorrectly marked as neuronal
FN = sum(true_label == 1 & pred_label == 0);  % Missed neuronal
TN = sum(true_label == 0 & pred_label == 0);  % Correctly marked non-neuronal

% Compute metrics
accuracy = (TP + TN) / length(true_label);
precision = TP / (TP + FP + eps);
recall = TP / (TP + FN + eps);
f1_score = 2 * (precision * recall) / (precision + recall + eps);



%confusion matrix
confusionMat = [TP, FN; FP, TN];
disp('Confusion Matrix:');
disp(array2table(confusionMat, ...
    'VariableNames', {'Predicted_Neuronal', 'Predicted_NonNeuronal'}, ...
    'RowNames', {'Actual_Neuron', 'Actual_NonNeuron'}));

%metrics
fprintf('\nMetrics:\n');
fprintf('True Positives (TP): %d\n', TP);
fprintf('False Positives (FP): %d\n', FP);
fprintf('False Negatives (FN): %d\n', FN);
fprintf('True Negatives (TN): %d\n', TN);
fprintf('Accuracy: %.4f\n', accuracy);
fprintf('Precision: %.4f\n', precision);
fprintf('Recall: %.4f\n', recall);
fprintf('F1 Score: %.4f\n', f1_score);

%metrics of match or no match to be imported in PHY
resultTable = group1(:, {'cluster_id'});  
resultTable.NvsNN_fp_fn = merged.fp_fn;
outFile = fullfile(pwd, sprintf('NvsNN_%s.tsv', username));
writetable(resultTable, outFile, 'Delimiter', '\t', 'FileType', 'text');
end