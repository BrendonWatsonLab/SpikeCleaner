%%%%%%%%%this code is for comparing good vs rest
 %user labels [rename cluster_group.tsv from phy to cluster_group.Username.tsv]
 %change its column header from 'group' to 'Usernamegroup'
group = readtable("cluster_group_Username.tsv", 'FileType', 'text', 'Delimiter', '\t');
group1 =readtable("path to cluster_TIDEV1.tsv", 'FileType', 'text', 'Delimiter', '\t');%algorithm

% Cleaning up strings
group.group = strtrim(lower(string(group.Usernamegroup))); %add the column header of user: from 'group-->groupUsername'
group1.reason = strtrim(lower(erase(erase(string(group1.TideV1), '('), ')')));

% Merge tables on cluster_id
merged = innerjoin(group, group1, 'Keys', 'cluster_id');


 %: "neuronal" if good or mua
true_label = ismember(merged.group, ["good", "mua"]);

% So, we mark "neuronal" as 1 if it's NOT labeled as noise/lowrate
non_neuronal_algo = contains(merged.reason, "noise") | contains(merged.reason, "lowrate");
pred_label = ~non_neuronal_algo;

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
    'VariableNames', {'Predicted Neuronal', 'Predicted Non-Neuronal'}, ...
    'RowNames', {'Actual Neuron', 'Actual Noise'}));

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