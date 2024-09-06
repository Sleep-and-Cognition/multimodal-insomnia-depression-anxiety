%% Create overview table with all brain-behavior associations

% Initialize an empty table
res_table = table();

% Assign the first column of the table with region descriptions
res_table(:,1) = regionDescriptions;

% Select symptom ('A', 'D', or 'I')
symptom = 'A';

% Define the data modalities
data_modalities = {'surfaces', 'thicknesses', 'volumes', 'connvS', 'connvF'};

% Iterate over each data modality
for id = 1:length(data_modalities)

    % Get the current data modality
    data_modality = data_modalities{id};

    % Construct the path to the results file
    results_file_path = fullfile(results_path_general, ...
        'main_analyses', [data_modality '*'], 'both');
    results_file_path = dir(results_file_path);
    results_file_path = results_file_path(1).folder;

    % Load the results from the file
    R = load(fullfile(results_file_path, 'results_script_present_regional_maps'));

    % Find the indices of the regions in the table
    [~, indx] = ismember(R.significance_regions.A(:,1), regionDescriptions);

    % Iterate over each region
    for i = 1:length(R.significance_regions.A(:,1))

        % report either the exact p-value or p < 0.001, or p = 1.000
        if R.significance_regions.(symptom){i,3} < 0.001
            res_table(indx(i), id+1) = {sprintf('β = %.3f, (p < 0.001)', ...
                R.significance_regions.(symptom){i,2})};
        elseif isnan(R.significance_regions.(symptom){i,3})
            res_table(indx(i), id+1) = {sprintf('β = %.3f, (p = 1.000)', ...
                R.significance_regions.(symptom){i,2})};
        else
            res_table(indx(i), id+1) = {sprintf('β = %.3f, (p = %.3f)', ...
                R.significance_regions.(symptom){i,2}, ...
                R.significance_regions.(symptom){i,3})};
        end
    end

end

% Display the resulting table
disp(res_table)
