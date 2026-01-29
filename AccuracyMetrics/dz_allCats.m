function dz_allCats(curdir)
%%By Diksha Zutshi
%run first: to chnage the shame of the data folder to make it compatible.
if nargin < 1 || isempty(curdir)
    curdir=pwd;
end 


%% This code will compare all categories: Noise,Good, MUA
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

% Assign "Match" or "NotMatch"
merged.(newVar) = string(merged.(newVar));
merged.SpikeCleaner = string(lower(merged.SpikeCleaner)); % predictions
% Combine lowrate -> noise
merged.SpikeCleaner(merged.SpikeCleaner == "lowrate") = "noise";

% Initialize fp_fn column
merged.fp_fn = repmat("NotMatch", height(merged), 1);

% Assign "Match" where equal
mask = merged.(newVar) == merged.SpikeCleaner;
merged.fp_fn(mask) = "Match";


% === CONFUSION MATRIX COMPONENTS
label_map = containers.Map( ...
    string({'good','mua','noise'}), ...
    [1,      2,     3] );

% === Allocate vectors ===
true_label = zeros(height(merged),1);
pred_label = zeros(height(merged),1);

% === Convert to numeric for confusion matrix ===
for i = 1:height(merged)
    true_label(i) = label_map(merged.(newVar)(i));
    pred_label(i) = label_map(merged.SpikeCleaner(i));
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
%metrics of match or no match to be imported in PHY
resultTable = group1(:, {'cluster_id'});  
resultTable.allcat_fp_fn = merged.fp_fn;


outFile = fullfile(pwd, sprintf('allcat_%s.tsv', username));
writetable(resultTable, outFile, 'Delimiter', '\t', 'FileType', 'text');

end

