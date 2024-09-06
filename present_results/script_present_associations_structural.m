%% Script to present the results of regional brain measures.

%% Initialization
cm = flipud(cbrewer('div', 'RdBu', 1000));

for imodality =  1:length(data_modalities)
    caxis_lim = caxisLimList(imodality);
    data_modality = data_modalities{imodality};

    fprintf('--- %s ---\n', data_modality);

    for istrat= 1:length(stratifications)

        % Setup all directories
        stratification = stratifications{istrat};

        results_path = fullfile(results_path_general, ...
            experimentName, data_modality, stratification);
        figures_path = fullfile(figures_path_general, ...
            experimentName, data_modality, stratification);

        checkdir(results_path);
        checkdir(figures_path);

        % Load the computed brain maps
        load(fullfile(results_path, ...
            ['results_all_', data_modality '_' stratification resType]));

        difference_t = resall; % the '_t' does not imply anything.
        difference_t_symptoms = resall;
        difference_t_perm = resall_perm;
        difference_t_bt = resall_bt;

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

        % Find significant regions
        significance_regions = [];
        significance_regions_uncorrected = [];
        significance_total = [];
        difference_p_bt = nan(nregions, nboot, nsymptoms);

        for isymptom = 1:nsymptoms

            difference_p = nan(nregions, 1);
            for ir = 1:nregions
                [~, difference_p(ir)] = ttest2(difference_t(ir, isymptom), ...
                    squeeze(difference_t_perm(ir, isymptom, :)));
            end

            % Calculate regional p-values in each bootstrap sample
            for ir = 1:nregions
                x = squeeze(difference_t_bt(ir, isymptom, :))';
                y = squeeze(difference_t_perm(ir, isymptom, :));
                n = nperm;
                dfe = n-1;
                s2y = nanvar(y);
                sPooled = sqrt(s2y);
                se = sPooled .* sqrt(1+ (1 ./ n));
                t = (mean(y,1) - x) ./ se;
                difference_p_bt(ir, :, isymptom) = 2 * tcdf(-abs(t),dfe);
            end

            % Calculate the FDR-corrected p-values (named 'a').
            [~, ~, a] = fdr(difference_p);

            % Present the signficant brain regions (sorted by effect size)
            regions_significantFDR = regionDescriptions(a < 0.05);
            [effects_significantFDR_sorted, I] = sort(...
                difference_t(a < 0.05, isymptom), 'descend');

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
            [~, total_difference_p] = ttest2(total_resall(isymptom), ...
                squeeze(total_resall_perm(1, isymptom, :)));
            fprintf('Total score: beta = %.3f, p = %.3g.\n', ...
                total_resall(:, isymptom), total_difference_p);

            % Save the (FDR-corrected) p-values such that they can later be
            % compared across experiments.
            significance_regions.(symptoms{isymptom}) = ...
                [regionDescriptions, ...
                num2cell(difference_t(:,isymptom)), ...
                num2cell(a)];

            significance_regions_uncorrected.(symptoms{isymptom}) = ...
                [regionDescriptions, ...
                num2cell(difference_t(:,isymptom)), ...
                num2cell(difference_p)];

            significance_total.(symptoms{isymptom}) = ...
                ['global', ...
                {total_resall(isymptom)}, ...
                {total_difference_p}];

            % Show the regional map on the brain
            plotBrain([regionDescriptions; {'extra'}], ...
                [difference_t(:,isymptom); 10], cm, ...
                'limits', [-0.05 .05], ...
                'savePath', fullfile(figures_path, symptoms{isymptom}), ...
                'atlas', [atlas '_aseg_column'], ...
                'Viewer', false);

            plotBrain([regionDescriptions; {'extra'}], ...
                [(a < 0.05) .* difference_t(:,isymptom); 10], cm, ...
                'limits', [-0.05 .05], ...
                'savePath', [fullfile(figures_path, ...
                symptoms{isymptom}) '_significant'], ...
                'atlas', [atlas '_aseg_column'], ...
                'Viewer', false);

        end

        close all

        save(fullfile(results_path, 'results_script_present_regional_maps'), ...
            'significance_regions', 'significance_total', ...
            'significance_regions_uncorrected', 'difference_p_bt', ...
            'difference_t_perm', 'difference_t_bt');

    end

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

    %% Follow-up analysis
    script_neurotransmission
    script_nke
    script_yeo
    
end