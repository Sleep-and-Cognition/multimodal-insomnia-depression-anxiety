%% Script to calculate the correlation of the association maps across modalities.
% NOTE: This script uses the variable spins, which is defined in
% script_use_spin_null.m

%% Initialize
delete(fullfile(results_path_general, 'crossmodal.txt'))
diary(fullfile(results_path_general, 'crossmodal.txt'));

load('regionDescriptionsFull.mat');

resTable = table();
res = zeros(82,3);
iExperiment = 1;

modality_names = {'surfaces', ...
    'thicknesses', ...
    'volumes', ...
    'reactivity', ...
    'SC', ...
    'FC'};

figure('color', 'white');
tiledlayout(1,3);

symptoms = {'I', 'D', 'A'};
symptomsFull = {'Insomnia symptoms', ...
    'Depressive symptoms', ...
    'Anxiety symptoms'};

for isymptom = 1:3
    %% Load regional results from all modalities
    symptom = symptoms{isymptom};
    symptomFull = symptomsFull{isymptom};
    symptomNum = isymptom;
    R = nan(length(regionDescriptions_all), 1); % temporary array

    % Load data
    resall = nan(size(regionDescriptions_all,1),6);
    pall = nan(size(regionDescriptions_all,1),6);

    results_path_experiment = fullfile(results_path_general, ...
        experimentalSettings(iExperiment).experimentName);

    results_path = fullfile(results_path_experiment, ['surfaces' '*'], stratification);
    results_path = dir(results_path); results_path = results_path(1).folder;
    R = load(fullfile(results_path, 'results_script_present_regional_maps'));
    resall(15:end, 1) = [R.significance_regions_uncorrected.(symptom){:,2}];
    pall(15:end, 1) = [R.significance_regions.(symptom){:,3}];
    resall(15:end, 1) = [R.significance_regions_uncorrected.(symptom){:,2}];

    results_path = fullfile(results_path_experiment, ['thicknesses' '*'], stratification);
    results_path = dir(results_path); results_path = results_path(1).folder;
    R = load(fullfile(results_path, 'results_script_present_regional_maps'));
    resall(15:end, 2) = [R.significance_regions_uncorrected.(symptom){:,2}];
    pall(15:end, 2) = [R.significance_regions.(symptom){:,3}];

    results_path = fullfile(results_path_experiment, ['volumes' '*'], stratification);
    results_path = dir(results_path); results_path = results_path(1).folder;
    R = load(fullfile(results_path, 'results_script_present_regional_maps'));
    resall(1:14,3) = [R.significance_regions_uncorrected.(symptom){:,2}];
    pall(1:14, 3) = [R.significance_regions.(symptom){:,3}];

    results_path = fullfile(results_path_experiment, ['connvS' '*'], stratification);
    results_path = dir(results_path); results_path = results_path(1).folder;
    R = load(fullfile(results_path, 'results_script_present_regional_maps'));
    resall(:,5) = [R.significance_regions_uncorrected.(symptom){:,2}];
    pall(:, 5) = [R.significance_regions.(symptom){:,3}];

    results_path = fullfile(results_path_experiment, ['connvF' '*'], stratification);
    results_path = dir(results_path); results_path = results_path(1).folder;
    R = load(fullfile(results_path, 'results_script_present_regional_maps'));
    resall(:,6) = [R.significance_regions_uncorrected.(symptom){:,2}];
    pall(:, 6) = [R.significance_regions.(symptom){:,3}];

    % tfmri
    results_path = fullfile(results_path_experiment, ['tfmri' '*'], stratification);
    results_path = dir(results_path); results_path = results_path(1).folder;
    R = load(fullfile(results_path, 'results_all_tfmri_both'));
    tfmri = R.resall;

    resall([6 13], 4) = tfmri(:,symptomNum);
    [~, pall([6 13], 4)] = ttest2(tfmri(symptomNum), squeeze(R.resall_perm(1,symptomNum,:)));

    R = arrayfun(@(x,p) sprintf('β = %.3f (p = %.3f)',x, p), resall, pall, 'UniformOutput', false);
    R = strrep(R, 'β = NaN (p = NaN)', 'NA');
    R = [regionDescriptionsFull R];

    %% Calculate correlations
    C = corr(resall, 'rows', 'pairwise');
        fprintf('Results for %s\n', symptomFull);
        disp(array2table(C, 'VariableNames',modality_names, ...
        'RowNames', modality_names));

    % Exclude t-fmri modality
    C = C([1 2 3 5 6], [1 2 3 5 6]);

    %% Use spin-based permutation testing to assess significance

    % Compute p-value using spin-based permutation testing
    resall_spin = nan(size(resall, 1), size(resall, 2), size(spins,2));

    % NOTE: the variable spins is defined in script_use_spin_null.m
    for ir = 1:size(resall,2)
        for i = 1:size(spins, 2)
            resall_spin(:, ir, i) = resall(spins(:, i), ir);
        end
    end

    C_spin = nan(size(resall,2), size(resall,2), size(spins,2));
    for i = 1:size(spins, 2)
        C_spin(:, :, i) = corr(resall, resall_spin(:, :, i), ...
            'rows', 'pairwise');
    end

    % Exclude t-fmri modality
    C_spin = C_spin([1 2 3 5 6], [1 2 3 5 6], :);

    p_spin = sum(abs(C_spin) >= abs(C), 3) ./ size(spins,2);
    p_spin(isnan(C)) = NaN;
    p_spin = setdiag(p_spin,0);
    p_spin = squareform(p_spin);

    % Apply FDR-correction
    [~, ~, a] = fdr(p_spin);
    pval_a = squareform(a);

    % Present significant correlations
    fprintf('FDR corrected across %i tests.\n', nnz(~isnan(p_spin)));    
    fprintf('p-values:\n');
    disp(array2table(pval_a, 'VariableNames',modality_names([1 2 3 5 6]), ...
        'RowNames', modality_names([1 2 3 5 6])));


    %% Present results in a figure

    % Prepare correlation matrix for presentation
    C(isnan(C)) = -1;
    C = setdiag(C,1);
    C(triu(C)~=0) = -1;

    a = nexttile;
  
    colormap(cm);
    clim([-1 1]);
    hold on;

    % Plot the correlation matrix "upside down" for visual reasons.
    C = flipud(C);
    pval_a = flipud(pval_a);

    % We cannot use the built-in imagesc function, because it does not
    % allow for hatching. We will use patch instead.
    % Compute the x and y coordinates of the patches
    nx = size(C, 1);
    ny = size(C,2);
    u = repmat([1:nx]', 1, ny);
    v = repmat(1:ny, nx, 1);
    
    clear h
    for i = 1:numel(C)
        y = u(i) + [-0.5 -0.5 0.5 0.5];
        R = v(i) + [-0.5 0.5 0.5 -0.5];
        h(i) = patch(R, y, C(i), ...
            'EdgeColor', [1 1 1], ...
            'LineWidth', 0.5);
    end

    for i = 1:length(h)
        if (pval_a(i) >= 0.05) || isnan(pval_a(i))
            hatch(h(i))
        end
    end
   
    title(symptomFull);
    axis square;

    a = gca;
    a.XTick = [];
    a.XAxis.Color = [ 1 1 1];
    a.YTick = [];
    a.YAxis.Color = [ 1 1 1];

    set(a, ...
        'FontName', 'Arial', ...
        'FontSize', 10, ...
        'Box', 'off', ...
        'LineWidth', 1);

    % Labels y-axis
    if tilenum(a) == 1
        str = modality_names([2 3 5 6]);
        str = fliplr(str); % because we put the matrix "upside down".
        text(zeros(1,4), 1:4, str, 'HorizontalAlignment', 'right');
    end

    % Labels x-axis
    y =  zeros(1,4);
    str = modality_names([1 2 3 5]);
    text(1:4, y, str, 'HorizontalAlignment', 'right', 'rotation', 45);

end

c = colorbar;

diary('off');

