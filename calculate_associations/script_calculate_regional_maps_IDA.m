%% Calculate symptom maps and null-model
% This script calculates the symptom maps and null-models for all
% modalities.

for imodality = 1:length(data_modalities)
    
    % Load the current data modality as `thisData`
    data_modality = data_modalities{imodality};
    thisData = eval(data_modality);
    disp(data_modality);
    
    % If the data modality is a surface or thickness, filter the data to
    % only include cortical regions. If the data modality is a volume,
    % include only subcortical regions. For other modalities (SC, FC,
    % tfmri) keep all regions.
    if any(contains(data_modality, {'surfaces', 'thicknesses'}))
        indx = contains(regionDescriptions_all, 'ctx');
        thisData = thisData(indx, :);
        regionDescriptions = regionDescriptions_all(indx);
    elseif any(contains(data_modality, {'volumes'}))
        indx = ~contains(regionDescriptions_all, 'ctx');
        thisData = thisData(indx, :);
        regionDescriptions = regionDescriptions_all(indx);
    else
        regionDescriptions = regionDescriptions_all;
    end

    if isempty(thisData)
        continue
    end
    
    % Add a row with the sum of all measures (total volume etc.).
    % This is equivalent to taking the average (therefore also used for
    % thicknesses).
    totalData = nansum(thisData,1);
    totalData(totalData == 0) = NaN;
    thisData = [thisData; totalData];
    
    % Get the number of regions
    nregions = length(regionDescriptions);
    
    for istrat= 1:length(stratifications)
        
        stratification = stratifications{istrat};
        
        % Create the results and figures paths for the current modality and
        % stratification
        results_path = fullfile(results_path_general, ...
            experimentName, data_modality, stratification);
        checkdir(results_path);
        
        figures_path = fullfile(figures_path_general, ...
            experimentName, data_modality, stratification);
        checkdir(figures_path);
        
        delete(fullfile(results_path, 'script_calculate_regional_maps_IDA.txt'));
        diary(fullfile(results_path, 'script_calculate_regional_maps_IDA.txt'));
        
        % Make symptom maps and permutations
        script_make_covariates
        
        % Exclude ovariates with no variance
        indx = nanstd(covariates(~all(isnan(thisData),1),:)) == 0;
        if any(indx)
            warning(['Problem with covariate numbers %s.\n', ...
                'Covariate has zero standard deviation and is not included as covariate.\n\n'], ...
                num2str(find(indx)));
        end
        covariates = covariates(:, ~indx);
        
        if 
        script_analysis_phenotype
        [resall, resall_perm, resall_bt] = getResults(thisData, YY, covariatesAll, nperm, nboot, nsample);
        
        % Alternative for debugging:
        % indx_others = ~any(isnan(YY),2) & ~any(isnan(covariates),2) & ~isnan(thisData(15,:))';
        % fitlm(zscore([thisData(1,indx_others)' covariates(indx_others,:)]), zscore(YY(indx_others,1)))
        
        % Extract results for total effects
        total_resall = resall(end,:);
        total_resall_perm = resall_perm(end,:,:);
        total_resall_bt = resall_bt(end,:,:);
        
        % Remove the last row from the result matrices (total effects)    
        resall = resall(1:end-1, :);
        resall_perm = resall_perm(1:end-1, :, :);
        resall_bt = resall_bt(1:end-1, :, :);
        
        % Save the results to a file
        save(fullfile(results_path, ['results_all_', data_modality '_' stratification resType]), ...
            'total_resall', 'total_resall_perm', 'resall', 'resall_perm', 'resall_bt', 'total_resall_bt', ...
            'regionDescriptions', 'nregions', 'symptoms', 'nsymptoms', 'indx_conn');
        
        diary('off');

    end
end