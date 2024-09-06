%% Add default covariates

fprintf('Default covariates:\n');

age = info_table.(['f_21003_' timePointImaging '_0']);
gender = info_table.f_31_0_0;

% Stratification option was not explored / used in this study.
switch stratification
    case 'both'
        indx_strat = true(size(gender));
        covariates = [age gender == 'Male'];
        fprintf('\tage\n\tgender\n');
end

covariates = [covariates age.^2 (age .* (gender == 'Male'))];
fprintf('\tage^2\n');
fprintf('\tage*gender\n');

covariates = [covariates covar_table.(['f_12144_' num2str(timePointImaging) '_0'])];
fprintf('\theight\n');

covariates = [covariates covar_table.('f_1707_0_0') == 'Right-handed'];
fprintf('\tRight handed\n');

%% Add imaging covariates

covariates_imaging_orig = covar_table(:, {['f_25756_' timePointImaging '_0'], ...
    ['f_25757_' timePointImaging '_0'], ['f_25758_' timePointImaging '_0'], ...
    ['f_25759_' timePointImaging '_0'], ['f_25926_' timePointImaging '_0']});

% Present histogram of items that are not excluded (i.e. NaNs)
indx = ~all(isnan(thisData),1);
figure('color', 'white');
for i = 1:size(covariates_imaging_orig,2)
    subplot(ceil(sqrt(size(covariates_imaging_orig,2))), ...
        round(sqrt(size(covariates_imaging_orig,2))), i);
    histogram(covariates_imaging_orig{indx, i});
    xlabel(textwrap(covariates_imaging_orig.Properties.VariableDescriptions(i), 75));
end

print(fullfile(figures_path, 'covariates_histogram'), '-dsvg');

fprintf('Imaging covariates:\n');
disp(covariates_imaging_orig.Properties.VariableDescriptions');

% Add assessment center
assessment_center = covar_table(:, ['f_54_' timePointImaging '_0']);
assessment_center_dummy = dummyvar(categorical(assessment_center{:,1}));
assessment_center_dummy(isnan(assessment_center{:,1})) = NaN;

% Exclude assessment centers with no participants
covariates_imaging = [covariates_imaging_orig{:, :} 
    assessment_center_dummy(:, 2:end)];

fprintf('\tAssessment centre:\n');
disp(grpstats(assessment_center, ['f_54_' timePointImaging '_0'], 'numel'));
