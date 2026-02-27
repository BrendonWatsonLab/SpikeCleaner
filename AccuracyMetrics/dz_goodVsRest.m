function dz_goodVsRest(curdir)
%%By Diksha Zutshi
%run first: to chnage the name of the data folder to make it compatible.
if nargin < 1 || isempty(curdir)
    curdir=pwd;
end    
%%%%%%%%%this code is for comparing good vs rest: SU Vs Rest
%%this will export a metrics.tsv for match or no match
username=strtrim(lower(input('Enter the user name for the PHY display:','s')));
group1 =readtable("cluster_SpikeCleaner.tsv", 'FileType', 'text', 'Delimiter', '\t');%spike cleaner outputs
userlabels='cluster_group.tsv';
userlabelpath=fullfile(curdir,userlabels);
renameduser=sprintf('cluster_group_%s.tsv',username);
renameduserpath=fullfile(curdir,renameduser);


if ~isfile(renameduserpath)
    group= readtable(userlabels, 'FileType', 'text', 'Delimiter', '\t');
    oldVar = 'group';
    newVar = sprintf('%s_group', username);
    if isfile(userlabelpath)
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
    
else
    group= readtable(renameduserpath, 'FileType', 'text', 'Delimiter', '\t');
end    


% Cleaning up strings
group1.SpikeCleaner = strtrim(lower(group1.SpikeCleaner));
% Merge tables on cluster_id
merged = innerjoin(group, group1, 'Keys', 'cluster_id');
%brendon_group
renamedcol=[username,'_','group'];



%labels: 1 = good, 0 = not good
true_label = double(merged.(renamedcol) == "good");
pred_label = double(merged.SpikeCleaner == "good");

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
resultTable.goodVsrest = repmat("", height(resultTable), 1);

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
%metrics of match or no match to be imported in PHY
outFile = fullfile(pwd, sprintf('GoodvsRest_%s.tsv', username));
writetable(resultTable, outFile, 'Delimiter', '\t', 'FileType', 'text');

end
