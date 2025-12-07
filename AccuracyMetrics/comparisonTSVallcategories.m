%% This code will compare all categories: Noise,Good, MUA

%user labels [rename cluster_group.tsv from phy to cluster_group.Username.tsv]
%change its column header from 'group' to 'Usernamegroup'
group = readtable("cluster_group_Username.tsv", 'FileType', 'text', 'Delimiter', '\t');
group1 =readtable("path to cluster_TIDEV1.tsv", 'FileType', 'text', 'Delimiter', '\t');%algorithm


% Clean up and normalize strings
group.group = strtrim(lower(string(group.Username)));
group1.reason = strtrim(lower(erase(erase(string(group1.TideV1), '('), ')')));
group1.reason(group1.reason == "lowrate") = "noise";  % unify lowrate into noise

% === MERGE TABLES BY CLUSTER ID ===
merged = innerjoin(group, group1, 'Keys', 'cluster_id');

% === ASSIGN MATCH ===
% Merge tables
merged = innerjoin(group, group1, 'Keys', 'cluster_id');

% Assign "Match" or "NotMatch"
merged.fp_fn = repmat("NotMatch", height(merged), 1);  % Default
merged.fp_fn(merged.group == merged.reason) = "Match"; % Assign match



% === CONFUSION MATRIX COMPONENTS (OPTIONAL)
label_map = containers.Map({'good', 'mua', 'noise'}, [1, 2, 3]);
true_label = zeros(height(merged), 1);
pred_label = zeros(height(merged), 1);

for i = 1:height(merged)
    true_label(i) = label_map(merged.group(i));
    pred_label(i) = label_map(merged.reason(i));
end

confusionMat = confusionmat(true_label, pred_label);
confusion_table = array2table(confusionMat, ...
    'VariableNames', {'Pred_Good', 'Pred_MUA', 'Pred_Noise'}, ...
    'RowNames', {'True_Good', 'True_MUA', 'True_Noise'});

disp('=== Confusion Matrix ===');
disp(confusion_table);

% === PRINT METRICS
TP = sum(true_label == 1 & pred_label == 1);
FP = sum(true_label ~= 1 & pred_label == 1);
FN = sum(true_label == 1 & pred_label ~= 1);
TN = sum(true_label ~= 1 & pred_label ~= 1);

accuracy = (TP + TN) / length(true_label);
precision = TP / (TP + FP + eps);
recall = TP / (TP + FN + eps);
f1_score = 2 * (precision * recall) / (precision + recall + eps);

fprintf('\n=== Metrics ===\n');
fprintf('Accuracy: %.4f\n', accuracy);
fprintf('Precision: %.4f\n', precision);
fprintf('Recall: %.4f\n', recall);
fprintf('F1 Score: %.4f\n', f1_score);
% % 
% % === EXPORT FINAL RESULT
outputFile ="path to Data folder /allcatUsername.tsv"; %metrics of match or no match
