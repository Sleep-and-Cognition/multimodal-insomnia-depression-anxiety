stratification = 'both';
data_modality = 'tfmri';
nsymptoms = 3;

%% Present tfmri
results_file_path = fullfile(results_path_general, ...
    experimentalSettings(iExperiment).experimentName, ...
    ['tfmri' '*'], stratification);

results_file_path = dir(results_file_path);
f = results_file_path(3).name;
results_file_path = results_file_path(1).folder;
data = load(fullfile(results_file_path, f));

figures_path = fullfile(figures_path_general, experimentName, ...
    data_modality, stratification);
results_path = fullfile(results_path_general, ...
    experimentName, data_modality, stratification);

p = nan(length(data.symptoms),1);
for is = 1:length(data.symptoms)
    [~,p(is)] = ttest2(data.resall(is), squeeze(data.resall_perm(1,is,:)));
end

[~, ~, a] = fdr(p);

% Setup logging
delete(fullfile(results_path, 'tfmri-general.txt'));
diary(fullfile(results_path, 'tfmri-general.txt'));

% Start reporting all settings
fprintf('experimentName:       %.80s\n', experimentName);
fprintf('timePointImaging:     %.80s\n', jsonencode(timePointImaging));
fprintf('timePoints:           %.80s\n', jsonencode(timePoints));
fprintf('covariatesAdditional: %.80s\n', jsonencode(covariatesAdditional));
fprintf('data_modality:        %.80s\n', data_modality);
fprintf('connectionsFlag:      %.80g\n', connectionsFlag);
fprintf('dichotomizeFlag:      %.80g\n', dichotomizeFlag);
fprintf('nperm:                %.80g\n', nperm);
fprintf('nboot:                %.80g\n\t', nboot);

fprintf('tfMRI results:\n');
fprintf('\tsymptom name\tEffect size\tFDR-corrected p-value\n');
res = [data.symptoms' num2cell(data.resall') num2cell(a)]';
fprintf('\t%s\t\t%.4f\t\t%.4f\n', res{:})



fprintf('Significant differences between global effects:\n');
res = nan(3);
for is = 1:nsymptoms
    for iis = is+1:nsymptoms
        x = squeeze(data.resall_bt(1,is,:) -  data.resall_bt(1,iis,:));
        [~,res(is,iis)] = ttest2(x,0);
    end
end

[~, ~, a] = fdr(res(:));
a = reshape(a, size(res));

for is = 1:nsymptoms
    for iis = is+1:nsymptoms
        if a(is, iis) < 0.05
            fprintf('\tSIGNIFICANT: %s - %s p = %.3f\n', symptoms{is}, symptoms{iis}, a(is,iis));
        else
            fprintf('\tNOT significant: %s - %s p = %.3f\n', symptoms{is}, symptoms{iis}, a(is,iis));
        end

    end
end

diary('off');

for isymptom = 1:nsymptoms
    plotBrain({'Left-Amygdala', 'extra'}, [data.resall(isymptom) 10], cm, ...
        'limits', [-0.05 .05], ...
        'savePath', fullfile(figures_path, ['tfmis-' symptoms{isymptom}]), ...
        'atlas', [atlas '_aseg'], ...
        'Viewer', true);
end