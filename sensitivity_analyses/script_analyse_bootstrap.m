%% Script for bootstrap validation analyses.
% The script computes the bootstrap samples and evaluates the correlation
% between the main analysis and alternative analyses with respect to the
% bootstrap samples.

%% Initialize variables
stratification = 'both';
symptoms = {'I', 'D', 'A'};
nboot = 5000;

%% Framework for running and presenting analyses across experimental settings

% Add specific text to the output file (not used here?)
resType = '';

stratifications = {'both'};

clear experimentalSettings

% Define the experimental settings corresponding to the main analyses
experimentalSettings(1).experimentName = 'main_analyses';
experimentalSettings(1).timePointImaging = '2'; 
experimentalSettings(1).timePoints = [1 3]; 
experimentalSettings(1).covariatesAdditional = [];
experimentalSettings(1).data_modalities = {'thicknesses', 'surfaces', 'volumes'};
experimentalSettings(1).connectionsFlag = false;
experimentalSettings(1).dichotomizeFlag = false;
experimentalSettings(1).indx_conn = [];
experimentalSettings(1).nperm = 5000;
experimentalSettings(1).nboot = 5000;
experimentalSettings(1).nullmodel = 'subjectLabelPermutation';

experimentalSettings(end+1).experimentName = 'main_analyses';
experimentalSettings(end).timePointImaging = '2';
experimentalSettings(end).timePoints = [1 3]; 
experimentalSettings(end).covariatesAdditional = [];
experimentalSettings(end).data_modalities = {'connvS'};
experimentalSettings(end).connectionsFlag = true;
experimentalSettings(end).dichotomizeFlag = false;
experimentalSettings(end).nperm = 5000;
experimentalSettings(end).indx_conn = indx_conn_SC;
experimentalSettings(end).nboot = 5000;
experimentalSettings(end).nullmodel = 'subjectLabelPermutation';

experimentalSettings(end+1).experimentName = 'main_analyses';
experimentalSettings(end).timePointImaging = '2'; 
experimentalSettings(end).timePoints = [1 3]; 
experimentalSettings(end).covariatesAdditional = [];
experimentalSettings(end).data_modalities = {'connvF'};
experimentalSettings(end).indx_conn = indx_conn_FC;
experimentalSettings(end).connectionsFlag = true;
experimentalSettings(end).nperm = 5000;
experimentalSettings(end).dichotomizeFlag = false;
experimentalSettings(end).nboot = 5000;
experimentalSettings(end).nullmodel = 'subjectLabelPermutation';

%% Compute bootstrap samples
% N: 1,976 (second imaging visit)
% N: 3,968  (Holdout)
% N: 15,357 (alternative symptom scores)
nsampleList = [15357]; %3968; 1976];

for iExperiment = 1:length(experimentalSettings)

    % Update the experimental settings
    experimentFields = fieldnames(experimentalSettings);
    for i = 1:length(experimentFields)
        eval([experimentFields{i} '= experimentalSettings(' ...
            num2str(iExperiment) ').' experimentFields{i} ';']);
    end

    disp(experimentalSettings(iExperiment));

    % Load the specific data modality, add experimentalSettings.name to
    % results and figures path and compute the connectivity and regional
    % associations for the different bootstrap samples. Run analyses
    nperm = 1;

    for isample = 1:length(nsampleList)
        nsample = nsampleList(isample);
        for imodality = 1:length(experimentalSettings(iExperiment).data_modalities)
            data_modalities = experimentalSettings(iExperiment).data_modalities(imodality);
            script_calculate_regional_maps_IDA

            save(fullfile(results_path_general, ...
                'main_analyses', data_modality, stratification, ...
                ['bootstrap_samples_' num2str(nsample)]), 'total_resall_bt', 'resall_bt');
        end
    end
    close all



end

%% Compute the regional results using the bootstrap samples
for isample = 1:length(nsampleList)

    nsample = nsampleList(isample);

    for iExperiment = 1:3

        % Update the experimental settings
        experimentFields = fieldnames(experimentalSettings);

        for i = 1:length(experimentFields)
            eval([experimentFields{i} '= experimentalSettings(' ...
                num2str(iExperiment) ').' experimentFields{i} ';']);
        end

        for imodality = 1:length(data_modalities)

            data_modality = data_modalities{imodality};
            disp(data_modality);
            stratification = stratifications{1};
            results_path = fullfile(results_path_general, ...
                experimentName, data_modality, stratification);

            % Load the computed brain maps
            load(fullfile(results_path, ...
                ['results_all_', data_modality '_' stratification resType]));
            load(fullfile(results_path, ['bootstrap_samples_' num2str(nsample)]));

            if ~connectionsFlag
                difference_t_pr_perm = resall_perm;
                difference_t_pr_bt = resall_bt;
            else
                difference_t_pr_bt = nan(size(indx_conn, 1), size(resall,2), nboot);
                difference_t_pr_bt(indx_conn,:,:) = resall_bt;

                difference_t_pr_perm = nan(size(indx_conn,1),size(resall,2),nperm);
                difference_t_pr_perm(indx_conn,:,:) = resall_perm;

                W = squareform_md(difference_t_pr_perm, 'tomatrix');
                difference_t_pr_perm = squeeze(nanmean(W,2));

                W = squareform_md(difference_t_pr_bt, 'tomatrix');
                difference_t_pr_bt = squeeze(nanmean(W,2));
            end

            % Compute the p-values for the regional results
            difference_p_bt = nan(nregions, nboot, nsymptoms);

            for isymptom = 1:nsymptoms

                % For each bootstrap sample we computed the p-values of the
                % regional brain - symptom associations. Each bootstrap
                % sample (1 value) is compared to the permuted samples
                % (5,000 values/permutations). Columns are the bootstrap
                % samples, rows are the regions, and the third dimension
                % are the symptoms (symptoms scores). We perform a two-sided
                % t-test to compare the bootstrap samples to the permuted
                % samples. The code is similar to ttest2 but implemented
                % here for computational reasons.

                for ir = 1:nregions
                    x = squeeze(difference_t_pr_bt(ir, isymptom, :))';
                    y = squeeze(difference_t_pr_perm(ir, isymptom, :));
                    dfe = length(y)-1; % nx + ny - 2 --> 1 + ny - 2
                    s2y = nanvar(x); 
                    sPooled = sqrt(s2y);
                    se = sPooled .* sqrt(1+ (1 ./ length(y)));
                    t = (mean(y,1) - x) ./ se;
                    difference_p_bt(ir, :, isymptom) = 2 * tcdf(-abs(t),dfe);
                end

                total_difference_p_bt = nan(nboot, 1);
                for i = 1:nboot
                    [~, total_difference_p_bt(i), stats] = ttest2(total_resall_bt(1,isymptom, i), ...
                        squeeze(total_resall_perm(1, isymptom, :)));
                end

                save(fullfile(results_path, ...
                    ['results_bootstrap_' num2str(nsample)]), ...
                    'difference_p_bt', 'difference_t_pr_perm', ...
                    'difference_t_pr_bt', 'total_difference_p_bt', ...
                    'total_resall_bt');

            end

        end

    end
end

%% Present regional validation and replication results

clear repSettings

% Select the settings applicable for the specific validation analysis
repSettings(1).data_modalities =  {'thicknesses', 'surfaces', 'volumes', 'connvF', 'connvS'};
repSettings(1).experiments = [1 2];
repSettings(1).nsample = 15357;
repSettings(1).rep_path = fullfile(absolute_path, 'results/240606/aparc');
experimentalSettings(1).experimentName = 'main_analyses';
experimentalSettings(2).experimentName = 'alt-pheno';

for irep = 1:length(repSettings)
    data_modalities = repSettings(irep).data_modalities;
    experiments = repSettings(irep).experiments;
    nsample = repSettings(irep).nsample;
    rep_path = repSettings(irep).rep_path;

    diary(fullfile(results_path_general, ...
        ['replication_' experimentalSettings(experiments(2)).experimentName '.txt']));

    symptoms_long = {'insomnia', 'depression', 'anxiety'};

    res_rep = [];
    res_norep = [];

    figure;

    for imodality = 1:length(data_modalities)
        data_modality = data_modalities{imodality};
        fprintf('---%s---\n', data_modality);

        iExperiment = experiments(1);
        results_file_path = fullfile(results_path_general, ...
            experimentalSettings(iExperiment).experimentName, ...
            [data_modality '*'], stratification);
        results_file_path = dir(results_file_path);
        results_file_path = results_file_path(1).folder;
        x = load(fullfile(results_file_path, 'results_script_present_regional_maps'));
        fprintf('%s vs ', experimentalSettings(iExperiment).experimentName);
        xx = load(fullfile(results_file_path, ['results_bootstrap_' num2str(nsample)]));

        iExperiment = experiments(2);
        results_file_path = fullfile(rep_path, ...
            experimentalSettings(iExperiment).experimentName, ...
            [data_modality '*'], stratification);
        results_file_path = dir(results_file_path);
        results_file_path = results_file_path(1).folder;
        xrep = load(fullfile(results_file_path, 'results_script_present_regional_maps'));
        fprintf('%s\n', experimentalSettings(iExperiment).experimentName);

        for isymptom = 1:length(symptoms)
            symptom = symptoms{isymptom};
            indx = [x.significance_regions.(symptom){:,3}] < 0.05;
            if nnz(indx) >0
                names = xrep.significance_regions_uncorrected.(symptom)(indx,1);
                p = [xrep.significance_regions_uncorrected.(symptom){indx,3}];
                b = [xrep.significance_regions_uncorrected.(symptom){indx,2}];
                [~, ~, a] = fdr(p);

                % Calculate reference
                a_bt = nan(length(a), nboot);
                for ib = 1:nboot
                    [~, ~, a_bt(:,ib)] = fdr(xx.difference_p_bt(indx, ib, isymptom));
                end

                a_pval = nnz(sum(a < 0.05) >= sum(a_bt < 0.05,1)) ./ nboot;
                a_expected = mean(sum(a_bt < 0.05,1));

                names = strrep(names, 'ctx-lh-', 'left ');
                names = strrep(names, 'ctx-rh-', 'right ');
                names = strrep(names, '-', ' ');

                r = ([names(a < 0.05) num2cell(b(a < 0.05)') num2cell(a(a<0.05)')])';
                fprintf(['%i out of %i regional effects replicated for ' ...
                    '%s symptoms (expected %.1f, range = [%.1f - %.1f], ' ...
                    'p = %.3f), '], ...
                    nnz(a < 0.05), nnz(indx), symptoms_long{isymptom}, ...
                    a_expected, prctile(sum(a_bt < 0.05,1), [5 100]), a_pval);

                v = ones(length(indx),1);
                v(indx) = b .* (a < 0.05);
                checkdir(fullfile(figures_path_general, 'sensitivity_analyses'));
                plotBrain([x.significance_regions.I(:,1); {'extra'}], ...
                    [v; 10], [cm; 0.7 0.7 0.7], ...
                    'atlas', [atlas '_aseg_column'], ...
                    'limits', [-0.05 .05], ...
                    'savePath', fullfile(figures_path_general, ...
                    'sensitivity_analyses', ...
                    ['replication_', ...
                    experimentalSettings(iExperiment).experimentName, '_', ...
                    data_modality, '_', ...
                    symptom]), ...
                    'Viewer', false);
            end

            subplot(length(data_modalities), 3, (imodality-1)*3+isymptom);
            hold on;

            r = corr([x.significance_regions_uncorrected.(symptom){:,2}]', ...
                squeeze(xx.difference_t_pr_perm(:,isymptom,:)), 'rows', 'pairwise');

            h = histogram(r, 'Normalization', 'probability', ...
                'FaceColor', [0.7 0.7 0.7], ...
                'EdgeColor', [1 1 1], ...
                'FaceAlpha', 0.5, ...
                'BinEdges', -1:0.1:1);

            r_bt = corr([x.significance_regions_uncorrected.(symptom){:,2}]', ...
                squeeze(xx.difference_t_pr_bt(:,isymptom,:)), 'rows', 'pairwise');

            histogram(r_bt, 'BinEdges', h.BinEdges, ...
                'Normalization', 'probability', ...
                'FaceColor', cmmetal(1,:), ...
                'EdgeColor', [1 1 1], ...
                'FaceAlpha', 0.5);

            r = corr([x.significance_regions_uncorrected.(symptom){:,2}]', ...
                [xrep.significance_regions_uncorrected.(symptom){:,2}]', 'rows', 'pairwise');

            line([r r], [0 max(ylim)], ...
                'color', cmmetal(1,:), ...
                'lineWidth', 1.5);

            xlabel('Correlation coefficient');
            ylabel('Frequency (normalized)');
            xlim([-1 1]);

            % p-value of the observed correlation (one-sided permutation testing)
            p = mean(r_bt <= r);

            fprintf(['observed r=%.2f, ' ...
                'expected r=%.2f [%.2f - %.2f] 95-CI\n'], ...
                r, mean(r_bt), prctile(r_bt, [5 100]));

            if p < 0.05
                stattext = sprintf(' p = %.3f', p);
            else
                stattext = '';
            end

            title([symptom ' - ' data_modality ' - ' ...
                experimentalSettings(iExperiment).experimentName stattext], ...
                'Interpreter', 'none');

            g = gcf;
            g.Units = 'centimeters';
            g.Position = [   35.2425    8.4667   48.4011   38.6997];

            print(fullfile(figures_path_general, ...
                ['replication_' experimentalSettings(iExperiment).experimentName]), ...
                '-dsvg');

        end

        set(gca, ...
            'FontName', 'Arial', ...
            'FontSize', 10, ...
            'Box', 'off', ...
            'LineWidth', 1);

        % Present replication global results:
        results_orig = [ x.significance_total.I{3} ...
            x.significance_total.D{3} ...
            x.significance_total.A{3} ];
        results_rep = [ xrep.significance_total.I{3} ...
            xrep.significance_total.D{3} ...
            xrep.significance_total.A{3} ];

        p = results_orig;
        indx = p < 0.05;
        if any(indx)
            p_rep = results_rep(indx);

            disp([symptoms(indx)' num2cell(p_rep)'])
        end

        % Present whether replication is in line with expectations.
        for i = 1:3
            if indx(i)
                % Updated expected range to make it a one-sided test.
                fprintf('observed beta=%.3f, expected beta=%.3f [%.3f - %.3f] 95-CI\n', ...
                    xrep.significance_total.(symptoms{i}){2}, ...
                    mean(squeeze(xx.total_resall_bt(1,i,:))), ...
                    prctile(squeeze(xx.total_resall_bt(1,i,:)), [5 100]));
            end
        end

        fprintf('Global efffects:\n');
        for symptom = symptoms

            if x.significance_total.(symptom{1}){3} < 0.05
                if xrep.significance_total.(symptom{1}){3} < 0.05
                    fprintf('Replicated: %s, beta=%.3f, p=%.3f\n', ...
                        symptom{1}, xrep.significance_total.(symptom{1}){2}, ...
                        xrep.significance_total.(symptom{1}){3});
                else
                    fprintf('NOT Replicated: %s, beta=%.3f, p=%.3f\n', ...
                        symptom{1}, xrep.significance_total.(symptom{1}){2}, ...
                        xrep.significance_total.(symptom{1}){3});
                end
            end
        end

        fprintf('\n');

    end

    diary('off');
end



