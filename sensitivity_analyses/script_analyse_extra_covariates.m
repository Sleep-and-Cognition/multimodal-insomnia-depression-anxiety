%% Analyze impact extra covariates on the main symptom outcomes

%% Initialize
bmi =  info_table.f_21001_2_0;
covariates_to_test = {'townsend', ...
    'education_years', ...
    'CTS', ...
    'ATS', ...
    'smoking', ...
    'diabetes', ...
    'hypertension', ...
    'alcohol', ...
    'bmi', ...
    'medication'};

stratification = 'both';
nsample = 0;

%% Compute the brain-behavior associations for each covariate

for iCovariates = 1:length(covariates_to_test)

    clear experimentalSettings
    this_covariate = covariates_to_test{iCovariates};
    experimentalSettings(1).experimentName = ['test_' this_covariate];
    experimentalSettings(1).timePointImaging = '2'; 
    experimentalSettings(1).timePoints = [1 3]; 
    experimentalSettings(1).covariatesAdditional = eval(this_covariate);
    experimentalSettings(1).data_modalities = {'thicknesses', 'surfaces', 'volumes'};
    experimentalSettings(1).connectionsFlag = false;
    experimentalSettings(1).dichotomizeFlag = false;
    experimentalSettings(1).nperm = 5000;
    experimentalSettings(1).nboot = 2;
    experimentalSettings(1).nullmodel = 'subjectLabelPermutation';
    experimentalSettings(1).caxisLimList = [0.04 0.04 0.04];

    experimentalSettings(end+1).experimentName = ['test_' this_covariate];
    experimentalSettings(end).timePointImaging = '2'; 
    experimentalSettings(end).timePoints = [1 3]; 
    experimentalSettings(end).covariatesAdditional = eval(this_covariate);
    experimentalSettings(end).data_modalities = {'connvS'};
    experimentalSettings(end).connectionsFlag = true;
    experimentalSettings(end).dichotomizeFlag = false;
    experimentalSettings(end).nperm = 5000;
    experimentalSettings(end).indx_conn = indx_conn_SC;
    experimentalSettings(end).nboot = 2;
    experimentalSettings(end).caxisLimList = [0.02];
    experimentalSettings(end).nullmodel = 'subjectLabelPermutation';


    experimentalSettings(end+1).experimentName = ['test_' this_covariate];
    experimentalSettings(end).timePointImaging = '2'; 
    experimentalSettings(end).timePoints = [1 3]; 
    experimentalSettings(end).covariatesAdditional = eval(this_covariate);
    experimentalSettings(end).data_modalities = {'connvF'};
    experimentalSettings(end).indx_conn = indx_conn_FC;
    experimentalSettings(end).connectionsFlag = true;
    experimentalSettings(end).nperm = 5000;
    experimentalSettings(end).dichotomizeFlag = false;
    experimentalSettings(end).nboot = 2;
    experimentalSettings(end).caxisLimList = [0.02];
    experimentalSettings(end).nullmodel = 'subjectLabelPermutation';

    %% Run experiments
    for iExperiment = 1:length(experimentalSettings)
        try

            % Update the experimental settings
            experimentFields = fieldnames(experimentalSettings);
            for i = 1:length(experimentFields)
                eval([experimentFields{i} '= experimentalSettings(' ...
                    num2str(iExperiment) ').' experimentFields{i} ';']);
            end

            disp(experimentalSettings(iExperiment));

            % Present demographics
            indxSubjects = ~all(isnan(volumes),1);
            n = nnz(indxSubjects);
            fprintf('Sample size: %i\n', n);

            age = info_table.(['f_21003_' timePointImaging '_0']);
            gender = info_table.f_31_0_0;
            r = mean(gender(indxSubjects) == 'Female');
            fprintf('Percentage female: %.0f%%\n', r*100);
            fprintf('Female N=%.0f, Male N=%.0f\n', ...
                nnz(gender(indxSubjects) == 'Female'), ...
                nnz(gender(indxSubjects) == 'Male'));

            fprintf('Median age = %i\n', median(age(indxSubjects)));
            fprintf('Age range = %i-%i\n', ...
                min(age(indxSubjects)), max(age(indxSubjects)));

            % Load the specific data modality

            % Run analyses (uncomment following line to run computations)
            % script_calculate_regional_maps_IDA

            close all

            if connectionsFlag
                data_modality = data_modalities{1};
                caxis_lim = caxisLimList;
                script_network_analyses
            elseif any(contains(data_modalities, 'tfmri'))
                script_present_tfmri
            else
                script_present_regional_maps
            end

        catch ME
            warning(jsonencode(ME));
        end

    end

end

%% Present correlation between the main analyses and  replication analyses
experimentalSettings(1).experimentName = 'main_analyses';
experimentalSettings(1).timePointImaging = '2'; 
experimentalSettings(1).timePoints = [1 3]; 
experimentalSettings(1).covariatesAdditional = [];
experimentalSettings(1).data_modalities = ...
    {'surfaces','thicknesses', 'volumes'};
experimentalSettings(1).connectionsFlag = false;
experimentalSettings(1).dichotomizeFlag = false;
experimentalSettings(1).indx_conn = [];
experimentalSettings(1).nperm = 5000;
experimentalSettings(1).nboot = 5000;
experimentalSettings(1).caxisLimList = [0.04 0.04 0.04];
experimentalSettings(1).nullmodel = 'subjectLabelPermutation';

clear res_overview
clear p_overview

for iCovariates = 1:length(covariates_to_test)

    this_covariate = covariates_to_test{iCovariates};
    experimentalSettings(2).experimentName = ['test_' this_covariate];
    experimentalSettings(2).timePointImaging = '2'; 
    experimentalSettings(2).timePoints = [1 3];
    experimentalSettings(2).covariatesAdditional = eval(this_covariate);
    experimentalSettings(2).data_modalities = ...
        {'thicknesses', 'surfaces', 'volumes', 'connvF', 'connvS'};
    experimentalSettings(2).connectionsFlag = false;
    experimentalSettings(2).dichotomizeFlag = false;
    experimentalSettings(2).nperm = 5000;
    experimentalSettings(2).nboot = 2;
    experimentalSettings(2).nullmodel = 'subjectLabelPermutation';
    experimentalSettings(2).caxisLimList = [0.04 0.04 0.04];

    % Update the experimental settings
    iExperiment = 1;
    experimentFields = fieldnames(experimentalSettings);
    for i = 1:length(experimentFields)
        eval([experimentFields{i} '= experimentalSettings(' ...
            num2str(iExperiment) ').' experimentFields{i} ';']);
    end

    disp(experimentalSettings(iExperiment));

    clear repSettings
    repSettings(1).data_modalities =  
    {'thicknesses', 'surfaces', 'volumes', 'connvF', 'connvS'};
    repSettings(1).experiments = [1 2];
    repSettings(1).rep_path = fullfile(absolute_path, 'results/240606/aparc');

    for irep = 1:length(repSettings)
        data_modalities = repSettings(irep).data_modalities;
        experiments = repSettings(irep).experiments;
        rep_path = repSettings(irep).rep_path;

        diary(fullfile(results_path_general, ['replication_' ...
            experimentalSettings(experiments(2)).experimentName '.txt']));

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

            % Load the computed brain maps
            load(fullfile(results_file_path, ...
                ['results_all_', data_modality '_' stratification resType]));

            if ~isfield(x, 'difference_t_pr_perm')
                x.difference_t_pr_perm = x.difference_t_perm;
                x.difference_t_pr_bt = x.difference_t_bt;
            end

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
                        [~, ~, a_bt(:,ib)] = fdr(x.difference_p_bt(indx, ib, isymptom));
                    end

                    a_pval = nnz(sum(a < 0.05) >= sum(a_bt < 0.05,1)) ./ nboot;
                    a_expected = mean(sum(a_bt < 0.05,1));

                    names = strrep(names, 'ctx-lh-', 'left ');
                    names = strrep(names, 'ctx-rh-', 'right ');
                    names = strrep(names, '-', ' ');

                    fprintf(['%i out of %i regional effects replicated ' ...
                        'for %s symptoms (expected %.1f, p = %.3f), '], ...
                        nnz(a < 0.05), nnz(indx), ...
                        symptoms_long{isymptom}, a_expected, a_pval);

                    res_overview.(data_modality){isymptom, iCovariates, 1} = nnz(a < 0.05);
                    res_overview.(data_modality){isymptom, iCovariates, 2} = nnz(indx);
                    p_overview.(data_modality){isymptom, iCovariates, 2} = a_pval;

                end

                subplot(length(data_modalities), 3, (imodality-1)*3+isymptom);
                hold on;

                r = corr([x.significance_regions_uncorrected.(symptom){:,2}]', ...
                    squeeze(x.difference_t_pr_perm(:,isymptom,:)), ...
                    'rows', 'pairwise');

                h = histogram(r, 'Normalization', 'probability', ...
                    'FaceColor', [0.7 0.7 0.7], ...
                    'EdgeColor', [1 1 1], ...
                    'FaceAlpha', 0.5, ...
                    'BinEdges', -1:0.1:1);

                r_bt = corr([x.significance_regions_uncorrected.(symptom){:,2}]', ...
                    squeeze(x.difference_t_pr_bt(:,isymptom,:)), ...
                    'rows', 'pairwise');

                histogram(r_bt, 'BinEdges', h.BinEdges, ...
                    'Normalization', 'probability', ...
                    'FaceColor', cmmetal(1,:), ...
                    'EdgeColor', [1 1 1], ...
                    'FaceAlpha', 0.5);

                r = corr([x.significance_regions_uncorrected.(symptom){:,2}]', ...
                    [xrep.significance_regions_uncorrected.(symptom){:,2}]', ...
                    'rows', 'pairwise');

                line([r r], [0 max(ylim)], ...
                    'color', cmmetal(1,:), ...
                    'lineWidth', 1.5);

                xlabel('Correlation coefficient');
                ylabel('Frequency (normalized)');
                xlim([-1 1]);

                p = mean(r_bt <= r);

                fprintf(['observed r=%.2f, expected r=%.2f ' ...
                    '[%.2f - %.2f] 95-CI\n'], r, mean(r_bt), ...
                    prctile(r_bt, [2.5 97.5]));

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

            end

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

            fprintf('Global efffects:\n');
            for symptom = symptoms

                if x.significance_total.(symptom{1}){3} < 0.05
                    if xrep.significance_total.(symptom{1}){3} < 0.05
                        fprintf('Replicated: %s, beta=%.3f, p=%.3f\n', ...
                            symptom{1}, ...
                            xrep.significance_total.(symptom{1}){2}, ...
                            xrep.significance_total.(symptom{1}){3});
                    else
                        fprintf('NOT Replicated: %s, beta=%.3f, p=%.3f\n', ...
                            symptom{1}, ...
                            xrep.significance_total.(symptom{1}){2}, ...
                            xrep.significance_total.(symptom{1}){3});
                    end
                end
            end



            fprintf('\n');


        end

    end

end

%% Show overview heatmap

res_variable = struct();
figure('color', 'white');
tiledlayout(5,1);
for data_modality = data_modalities
    nexttile
    Obs = cell2mat(res_overview.(data_modality{1})(:, :, 1));
    Ref = cell2mat(res_overview.(data_modality{1})(:, :, 2));
    R = Obs ./ Ref;
    P = cell2mat(p_overview.(data_modality{1})(:, :));
    
    % Add exception when not all possible associations were present in main
    % analyses
    if size(R, 1) == 2
        R = [R; zeros(1, size(R,2))];
        P = [P; ones(1, size(P,2))];
        Obs = [Obs; zeros(1, size(P,2))];
        Ref = [Ref; zeros(1, size(P,2))];
    end

    for i = 1:size(P,1)
        for ii = 1:size(P,2)
            res_variable.(data_modality{1}){i,ii} = sprintf('%i/%i (p=%.3f)', ...
                Obs(i,ii), ...
                Ref(i,ii), ...
                P(i,ii));
        end
    end

    res_variable_table.(data_modality{1}) = 
    cell2table(res_variable.(data_modality{1}), 'VariableNames', ...
        {'SES', ...
        'Educational attainment', ...
        'Childhood stress', ...
        'Adulthood stress', ...
        'Smoking', ...
        'Diabetes', ...
        'Hypertension', ...
        'Alcohol use', ...
        'BMI', ...
        'Medication use'}, ...
        'RowNames', ...
        {'Insomnia symptoms', ...
        'Depression symptoms', ...
        'Anxiety symptoms'});

    a = gca;

    h = imagesc(R);

    hold on;
    % Get indices where elements of P are less than 0.05
    [rows, cols] = find(P < 0.05);

    % Overlay stars on these indices
    plot(cols, rows, 'p', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');

    cm = (cbrewer('seq', 'Blues', 1000));
    cm(1,:) = [1 1 1];
    caxis([0 1]);
    colormap(cm);

    xticks(1:13);
    yticks(1:3);

    yticklabels({'Insomnia symptoms', 'Depression symptoms', 'Anxiety symptoms'});

    xticklabels({'SES', 'Educational attainment', ...
        'Childhood stress', 'Adulthood stress', 'Smoking', ...
        'Diabetes', 'Hypertension', 'Alcohol use', 'BMI', 'Medication use'});

    set(a, ...
        'FontName', 'Arial', ...
        'FontSize', 10, ...
        'Box', 'off', ...
        'LineWidth', 1);

    colorbar

    axis equal
    axis tight

end
