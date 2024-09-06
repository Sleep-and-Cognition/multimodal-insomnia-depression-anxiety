%% Present structural and functional connectivity results

%% Initialization
stratification = 'both';

results_path = fullfile(results_path_general, experimentName, ...
    data_modality, stratification);
figures_path = fullfile(figures_path_general, experimentName, ...
    data_modality, stratification);

checkdir(results_path);
checkdir(figures_path);

%% Load results
% indx_conn refers to either SC or FC indx and is defined in
% script_main_run_experiments.m

load(fullfile(results_path, ...
    ['results_all_', data_modality '_' stratification resType]));

difference_t = nan(size(indx_conn, 1), size(resall,2));
difference_t(indx_conn,:) = resall;

difference_t_bt = nan(size(indx_conn, 1), size(resall,2), nboot);
difference_t_bt(indx_conn,:,:) = resall_bt;

difference_t_perm = nan(size(indx_conn,1),size(resall,2),nperm);
difference_t_perm(indx_conn,:,:) = resall_perm;

% Calculate the average effect across connections per region
W = squareform_md(difference_t, 'tomatrix');
difference_t_pr = squeeze(nanmean(W, 2));

W = squareform_md(difference_t_perm, 'tomatrix');
difference_t_pr_perm = squeeze(nanmean(W,2));

W = squareform_md(difference_t_bt, 'tomatrix');
difference_t_pr_bt = squeeze(nanmean(W,2));

% Save the results in variables that will be stored on disk, and will be
% used in the analyses across experiments.
significance_regions = [];
significance_regions_uncorrected = [];
significance_total = [];
difference_p_bt = nan(nregions, nboot, nsymptoms);

switch nullmodel
    case 'subjectLabelPermutation'
        fprintf('Subject label permutation used.\n');
    case 'spintesting'
        fprintf('Spin testing permutation used.\n')
        script_use_spin_null

        results_path = fullfile(results_path, 'spin');
        checkdir(results_path);
        figures_path = fullfile(figures_path, 'spin');
        checkdir(figures_path);
    otherwise
        fprintf(['nullmodel parameter not found. ' ...
            'Subject label permutation assumed.\n']);
end

% Setup logging
delete(fullfile(results_path, 'general.txt'));
diary(fullfile(results_path, 'general.txt'));

% Start reporting all settings
fprintf('experimentName:       %.80s\n', experimentName);
fprintf('timePointImaging:     %.80s\n', jsonencode(timePointImaging));
fprintf('timePoints:           %.80s\n', jsonencode(timePoints));
fprintf('covariatesAdditional: %.80s\n', jsonencode(covariatesAdditional));
fprintf('data_modality:        %.80s\n', data_modality);
fprintf('connectionsFlag:      %.80g\n', connectionsFlag);
fprintf('dichotomizeFlag:      %.80g\n', dichotomizeFlag);
fprintf('nperm:                %.80g\n', nperm);
fprintf('nboot:                %.80g\n', nboot);

for isymptom = 1:nsymptoms

    % Calculate regional p-values
    difference_p = nan(nregions, 1);
    for ir = 1:nregions
        [~, difference_p(ir)] = ttest2(difference_t_pr(ir, isymptom), ...
            squeeze(difference_t_pr_perm(ir, isymptom, :)));
    end
    [~, ~, a] = fdr(difference_p);

    % Calculate regional p-values in each bootstrap sample
    for ir = 1:nregions
        x = squeeze(difference_t_pr_bt(ir, isymptom, :))';
        y = squeeze(difference_t_pr_perm(ir, isymptom, :));
        dfe = n-1;
        s2y = nanvar(y);
        sPooled = sqrt(s2y);
        se = sPooled .* sqrt(1+ (1 ./ n));
        t = (mean(y,1) - x) ./ se;
        difference_p_bt(ir, :, isymptom) = 2 * tcdf(-abs(t),dfe);
    end

    % Calculate p-value of total brain differences
    [~, total_difference_p] = ttest2(total_resall(isymptom), ...
        squeeze(total_resall_perm(1, isymptom, :)));

    % Present the signficant brain regions (sorted by effect size)
    regions_significantFDR = regionDescriptions(a < 0.05);
    [effects_significantFDR_sorted, I] = sort(...
        difference_t_pr(a < 0.05, isymptom), 'descend');

    a_significantFDR = a(a < 0.05);
    result = [regions_significantFDR(I) ...
        num2cell(effects_significantFDR_sorted) ...
        num2cell(a_significantFDR(I))]';

    fprintf('%s:\n',symptoms{isymptom});
    fprintf('regions significant:\n');
    fprintf('%-30.30s\tbeta = %.3f, p = %.3f\n', result{:})

    fprintf('total number of significant regions: %i\n', ...
        length(a_significantFDR));

    % Present significance of total brain effects:
    fprintf('Total score: beta = %.2g, p = %.2g.\n', ...
        total_resall(:, isymptom), total_difference_p);

    % Save the (FDR-corrected) p-values such that they can later be
    % compared across experiments.
    significance_regions.(symptoms{isymptom}) = ...
        [regionDescriptions, ...
        num2cell(difference_t_pr(:,isymptom)), ...
        num2cell(a)];

    significance_regions_uncorrected.(symptoms{isymptom}) = ...
        [regionDescriptions, ...
        num2cell(difference_t_pr(:,isymptom)), ...
        num2cell(difference_p)];

    significance_total.(symptoms{isymptom}) = ...
        ['global', ...
        {total_resall(isymptom)}, ...
        {total_difference_p}];

    % Show the regional map on the brain
    % The 'extra' term added to make sure that all brain regions are white
    % in case they have all NaN values.
    plotBrain([regionDescriptions; {'extra'}], ...
        [difference_t_pr(:,isymptom); 10], cm, ...
        'savePath', fullfile(figures_path, [symptoms{isymptom}]), ...
        'limits', [-0.03 0.03], ...
        'atlas', [atlas '_aseg_column'], ...
        'Viewer', false);

    plotBrain([regionDescriptions; {'extra'}], ...
        [(a < 0.05) .* difference_t_pr(:,isymptom); 10], cm, ...
        'savePath', fullfile(figures_path, ['significant_' symptoms{isymptom}]), ...
        'limits', [-0.03 0.03], ...
        'atlas', [atlas '_aseg_column'], ...
        'Viewer', false);

end

save(fullfile(results_path, 'results_script_present_regional_maps'), ...
    'significance_regions', 'significance_total', ...
    'significance_regions_uncorrected', 'difference_p_bt', ...
    'difference_t_pr_bt', 'difference_t_pr_perm');

%% Compare global effects across disorders

fprintf('Significant differences between global effects:\n');
res = nan(3);
for is = 1:nsymptoms
    for iis = is+1:nsymptoms
        x = squeeze(total_resall_bt(1,is,:) -  total_resall_bt(1,iis,:));
        [~,res(is,iis)] = ttest2(x,0);
    end
end

[~, ~, a] = fdr(res(:));
a = reshape(a, size(res));

for is = 1:nsymptoms
    for iis = is+1:nsymptoms
        if a(is, iis) < 0.05
            fprintf('\tSIGNIFICANT: %s - %s p = %.3f\n', ...
                symptoms{is}, symptoms{iis}, a(is,iis));
        else
            fprintf('\tNOT significant: %s - %s p = %.3f\n', ...
                symptoms{is}, symptoms{iis}, a(is,iis));
        end

    end
end


diary('off');

%% Present network overviews

figure('Color', 'white');
for isymptom = 1:nsymptoms

    % Outcome network map
    W = squareform_md(difference_t, 'tomatrix');

    % Find the coordinates of the brain regions
    [~, I, J] = intersect(regionDescriptions, regionDescriptionsAll);
    W = W(I, I, isymptom);
    W = setdiag(W,nan);
    coordinates = regionCoordinates(J,:);

    W(abs(W(:)) < prctile(abs(W(~isnan(W))), 90)) = 0;
    W(isnan(W)) = 0;

    % Calculate the edge color based on their strength and colormap CM
    edge_color = squareform(W(:,:));
    edge_color = edge_color(edge_color ~= 0);
    edge_color = value2Color(edge_color, cm, -0.03, 0.03);

    % Represent nodes and edges as a graph. Matlab can plot graphs using "plot"
    G = graph(W(:, :) ~= 0);

    subplot(1,3,isymptom);
    h = plot(G, ...
        'XData', coordinates(:,1), ...
        'YData', coordinates(:,2), ...
        'EdgeColor', edge_color , ...
        'NodeColor', [0.7 0.7 0.7], ...
        'MarkerSize', 5, ...
        'LineWidth', 2, ...
        'NodeLabel', {});
    axis equal
    axis off

    title(symptoms(isymptom))

    colormap(cm);
    caxis([-0.3 0.3]);

    xlim([-70 70]);
    ylim([-100 100]);
    yticks([]);
    xticks([]);
    ax = gca;
    ax.Position = [(isymptom-1)*0.3 0.05 0.3 0.9];

end

g = gcf;
g.Units = 'normalized';

print(fullfile(figures_path, 'symptom_networks'), '-dsvg');

%% Follow-up analysis
script_neurotransmission
script_nke
script_yeo
