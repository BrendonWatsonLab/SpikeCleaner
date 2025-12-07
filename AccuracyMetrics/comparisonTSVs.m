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

%labels: 1 = good, 0 = not good
true_label = double(merged.group == "good");
pred_label = double(merged.reason == "good");

% Confusion matrix components
TP = sum(true_label == 1 & pred_label == 1);
FP = sum(true_label == 0 & pred_label == 1);
FN = sum(true_label == 1 & pred_label == 0);
TN = sum(true_label == 0 & pred_label == 0);

% Metrics
accuracy = (TP + TN) / length(true_label);
precision = TP / (TP + FP + eps);
recall = TP / (TP + FN + eps);
f1_score = 2 * (precision * recall) / (precision + recall + eps);

%confusion matrix
confusionMat = [TP, FN; FP, TN];
disp('Confusion Matrix:');
disp(array2table(confusionMat, ...
    'VariableNames', {'Predicted_Good', 'Predicted_NotGood'}, ...
    'RowNames', {'Actual_Good', 'Actual_NotGood'}));

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


%%%%%%%%%creating a new file according to cluster ids for comparison btw algo's and expert's decision
merged.fp_fn = repmat("", height(merged), 1);
merged.fp_fn(true_label == 1 & pred_label == 1) = "Match";  % TP
merged.fp_fn(true_label == 0 & pred_label == 0) = "Match";  % TN
merged.fp_fn(true_label == 0 & pred_label == 1) = "FP";  % FP
merged.fp_fn(true_label == 1 & pred_label == 0) = "FN";  % FN;

resultTable = group1(:, {'cluster_id'});  
resultTable.fp_fn = repmat("", height(resultTable), 1);

% index matchcing resultTable based on merged table  
for i = 1:height(merged)
    clusterID = merged.cluster_id(i);
    status = merged.fp_fn(i);
    
    % Find the matching row in the full resultTable
    rowIndex = find(resultTable.cluster_id == clusterID);
    if ~isempty(rowIndex)
        resultTable.fp_fn(rowIndex) = status;
    end
end
% % === EXPORT FINAL RESULT
writetable(resultTable, "Place  it in the data folder \metricsUsername.tsv", 'Delimiter', '\t', 'FileType', 'text');%metrics of match or no match
