%% Define covariates and symptom scores

covariatesAll = [ covariatesAdditional covariates];
covariatesAll = normalize(covariatesAll,1);

% Symptom severity scores

if strcmp(phenotype, 'alt-pheno')
    
    % Use alternative symptom severity scores

    % Insomnia symptom severity score
    Y1 = [nanmean(insomnia_score(:,timePoints),2)];

    % Depressive symptom severity score
    x = [phq2(:, timePoints) phq9_score depression_diagnosis ];
    x = normalize(x, "zscore");
    [~, SCORE] = pca(x, 'rows', 'complete');
    Y2 = SCORE(:,1);

    % Anxiety symptom severity score
    x = [n12_nervous(:, timePoints) GAD7 anx_diagnosis];
    x = normalize(x, "zscore");
    [~, SCORE] = pca(x, 'rows', 'complete');
    Y3 = SCORE(:,1);

else

    % Use default symptom severity scores
    
    Y1 = [nanmean(insomnia_score(:,timePoints),2)];
    Y2 = nanmean(phq2(:, timePoints),2);
    Y3 = nanmean(n12_nervous(:, timePoints),2);

end

if dichotomizeFlag
    Y1(~isnan(Y1)) = Y1(~isnan(Y1)) > 0;
    Y2(~isnan(Y2)) = Y2(~isnan(Y2)) > 0;
    Y3(~isnan(Y3)) = Y3(~isnan(Y3)) > 0;
end

YY = [Y1 Y2 Y3];
symptoms = {'I', 'D', 'A'};

nsymptoms = length(symptoms);

%% Report correlations between symptom severity scores
indxSubjects = ~all(isnan(thisData),1);
[r12, p12] = corr(Y1(indxSubjects), Y2(indxSubjects), ...
    "Rows","pairwise", "type","Spearman");
[r13, p13] = corr(Y1(indxSubjects), Y3(indxSubjects), ...
    "Rows","pairwise", "type","Spearman");
[r23, p23] = corr(Y3(indxSubjects), Y2(indxSubjects), ...
    "Rows","pairwise", "type","Spearman");

[~, ~, a] = fdr([p12 p23 p13]);

fprintf(['Spearman rank correlation coefficients: ' ...
    '%s-%s rho = %.2f, p = %.3f\n'], ...
    symptoms{1}, symptoms{2}, r12, a(1));
fprintf(['Spearman rank correlation coefficients: ' ...
    '%s-%s rho = %.2f, p = %.3f\n'], ...
    symptoms{2}, symptoms{3}, r23, a(3));
fprintf(['Spearman rank correlation coefficients: ' ...
    '%s-%s rho = %.2f, p = %.3f\n'], ...
    symptoms{1}, symptoms{3}, r13, a(2));