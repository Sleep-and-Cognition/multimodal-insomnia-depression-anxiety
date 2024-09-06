%% Present replication significant regions
diary(fullfile(results_path_general, 'sensitivity_excrelatives.txt'));
clear repSettings

repSettings(1).data_modalities =  {'thicknesses', 'surfaces', 'volumes'};
repSettings(2).data_modalities = {'connvF'};
repSettings(3).data_modalities = {'connvS'};

fprintf('----Family relations:----\n')

for irep = 1:length(repSettings)
data_modalities = repSettings(irep).data_modalities;

symptoms_long = {'insomnia', 'depression', 'anxiety'};

res_rep = [];
res_norep = [];

for imodality = 1:length(data_modalities)
data_modality = data_modalities{imodality};
fprintf('---%s---\n', data_modality);

results_file_path = fullfile(results_path_general, ...
    'main_analyses', ...
    [data_modality '*'], stratification);
results_file_path = dir(results_file_path);
results_file_path = results_file_path(1).folder;
x = load(fullfile(results_file_path, ...
    'results_script_present_regional_maps'));
fprintf('%s vs ', 'main_analyses');

results_file_path = fullfile(results_path_general, ...
    'excrelatives', ...
    [data_modality '*'], stratification);
results_file_path = dir(results_file_path);
results_file_path = results_file_path(1).folder;
xrep = load(fullfile(results_file_path, ...
    'results_script_present_regional_maps'));
fprintf('%s\n', 'excrelatives');


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
        
        r = ([names(a < 0.05) num2cell(b(a < 0.05)') ...
            num2cell(a(a<0.05)')])';

        fprintf(['%i out of %i regional effects replicated ' ...
            'for %s symptoms (expected %.1f, p = %.3f), '], ...
            nnz(a < 0.05), nnz(indx), ...
            symptoms_long{isymptom}, a_expected, a_pval);
       
    end
    
end
fprintf('\n');

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

diary('off');