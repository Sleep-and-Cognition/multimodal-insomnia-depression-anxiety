%% Validate the symptom severity scores

%% Initialization

% Compute depression and anxiety diagnostic data
diagnoses = {'Prefer not to answer (group B)', ...
        'Prefer not to answer (group A)', ...
        'Social anxiety or social phobia', ...
        'Schizophrenia', ...
        'Any other type of psychosis or psychotic illness', ...
        'A personality disorder', ...
        'Any other phobia (eg disabling fear of heights or spiders)', ...
        'Panic attacks', ...
        'Obsessive compulsive disorder (OCD)', ...
        'Mania, hypomania, bipolar or manic-depression', ...
        'Depression', ...
        'Bulimia nervosa', ...
        'Psychological over-eating or binge-eating', ...
        'Autism, Asperger''s or autistic spectrum disorder', ...
        'Anxiety, nerves or generalized anxiety disorder', ...
        'Anorexia nervosa', ...
        'Agoraphobia', ...
        ['Attention deficit or attention deficit and ' ...
        'hyperactivity disorder (ADD/ADHD)']};

for i = 1:16
    thisColumn = info_table.f_20544_0_1;
    iname = strcat('f_20544_0_', num2str(i));
    iname2 = strcat(iname, '_cat');
    info_table(:, iname2) = table(categorical(thisColumn, ...
        [-819, -818, 1, 2, 3, 4, 5, 6, 7, 10, 11, 12, 13, 14, 15, 16, 17, 18], ...
        {'Prefer not to answer (group B)', ...
        'Prefer not to answer (group A)', ...
        'Social anxiety or social phobia', ...
        'Schizophrenia', ...
        'Any other type of psychosis or psychotic illness', ...
        'A personality disorder', ...
        'Any other phobia (eg disabling fear of heights or spiders)', ...
        'Panic attacks', ...
        'Obsessive compulsive disorder (OCD)', ...
        'Mania, hypomania, bipolar or manic-depression', ...
        'Depression', ...
        'Bulimia nervosa', ...
        'Psychological over-eating or binge-eating', ...
        'Autism, Asperger''s or autistic spectrum disorder', ...
        'Anxiety, nerves or generalized anxiety disorder', ...
        'Anorexia nervosa', ...
        'Agoraphobia', ...
        ['Attention deficit or attention deficit and ' ...
        'hyperactivity disorder (ADD/ADHD)']}));
end

diagnose_fields = strcat('f_20544_0_', cellfun(@num2str, num2cell(1:16), ...
    'UniformOutput', 0));
diagnose_cat_fields = strcat('f_20544_0_', cellfun(@num2str, num2cell(1:16), ...
    'UniformOutput', 0), '_cat');

depression_diagnosis = double(any(info_table{:, ...
    diagnose_cat_fields} == 'Depression',2));
anx_diagnosis = double(any(info_table{:, ...
    diagnose_cat_fields} == ...
    'Anxiety, nerves or generalized anxiety disorder',2));

missing_items = any(info_table{:, diagnose_cat_fields} == ...
    'Prefer not to answer (group B)',2);
missing_items = missing_items | any(info_table{:, diagnose_cat_fields} == ...
    'Prefer not to answer (group A)',2);
missing_items = missing_items | ...
    any([info_table.f_20499_0_0 < 0 info_table.f_20500_0_0 < 0], 2);
missing_items = missing_items | any([ismissing(info_table.f_20499_0_0) ...
    ismissing(info_table.f_20500_0_0)], 2);

depression_diagnosis(missing_items) = NaN;
anx_diagnosis(missing_items) = NaN;

% Compute GAD-7 scores:
GAD7_items = [info_table.f_20505_0_0 ...
    info_table.f_20506_0_0 ...
    info_table.f_20509_0_0 ...
    info_table.f_20512_0_0 ...
    info_table.f_20515_0_0 ...
    info_table.f_20516_0_0 ...
    info_table.f_20520_0_0];

GAD7 = sum(GAD7_items - 1, 2);
GAD7_missing = any(GAD7_items < 0,2) | any(isnan(GAD7_items), 2);
GAD7(GAD7_missing) = NaN;

% Compute PHQ-9 and PHQ-2 scores:
depression_items = [info_table.f_20518_0_0 ...
    info_table.f_20507_0_0 ...
    info_table.f_20519_0_0 ...
    info_table.f_20511_0_0 ...
    info_table.f_20513_0_0 ...
    info_table.f_20508_0_0 ...
    info_table.f_20510_0_0 ...
    info_table.f_20514_0_0 ...
    info_table.f_20517_0_0];

phq9_score = sum(depression_items - 1, 2);
phq9_score(phq9_score < 0) = NaN;

phq2_score = sum(depression_items(:, [7 8 ]) - 1, 2);
phq2_score(phq2_score < 0) = NaN;

%% Calculate validity of PHQ2 with respect to depression diagnosis
fprintf('Depression validation:\n');

y = depression_diagnosis;
yhat = double(Y2 >= 1);
yhat(isnan(Y2)) = NaN;

TP = nnz((y == 1) & (yhat == 1));
TN = nnz((y == 0) & (yhat == 0));

FP = nnz((y == 0) & (yhat == 1));
FN = nnz((y == 1) & (yhat == 0));

% Sensitivity:
% TP / (TP + FN)
fprintf('Sensitivity: %.2f\n', TP / (TP + FN));

% Specificity
% TN / (FP + TN)
fprintf('Specificity: %.2f\n', TN / (FP + TN));

fprintf('Depression diag: N=%i, %.2f%%\n', ...
    nnz(~isnan(depression_diagnosis)), ...
    100*nnz(~isnan(depression_diagnosis)) ./ length(depression_diagnosis));


%% Validate PHQ2 versus PHQ9
[r,p] = corr(phq2_score, phq9_score, ...
    'rows', 'pairwise', ...
    'type', 'spearman');

fprintf('PHQ-9: N=%i, %.2f%%\n', ...
    nnz(~isnan(phq9_score)), ...
    100*nnz(~isnan(phq9_score)) ./ length(phq9_score));

%% Validate anxiety wrt. GAD7 or diagnosis
fprintf('Anxiety validation:\n');

% Validate with respect to GAD7:
% y = GAD7 > 15;

% Validate with respect to diagnosis:
y = anx_diagnosis;

yhat = double(Y3 > 0);
yhat(isnan(Y3)) = NaN;

TP = nnz((y == 1) & (yhat == 1));
TN = nnz((y == 0) & (yhat == 0));

FP = nnz((y == 0) & (yhat == 1));
FN = nnz((y == 1) & (yhat == 0));

% Sensitivity:
% TP / (TP + FN)
fprintf('Sensitivity: %.2f\n', TP / (TP + FN));

% Specificity
% TN / (FP + TN)
fprintf('Specificity: %.2f\n', TN / (FP + TN));



