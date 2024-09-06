%% Framework for running and presenting analyses across experimental settings

% Add specific text to the output file. 
resType = ''; % Not sure how this is used.

stratifications = {'both'}; % NOT USED

nsample = 0;

excludeRelatives = false;

clear experimentalSettings

% Main analysis
experimentalSettings(1).experimentName = 'main_analyses';
experimentalSettings(1).timePointImaging = '2'; % Imaging related data
experimentalSettings(1).timePoints = [1 3]; % Averaging over behavior related data
experimentalSettings(1).covariatesAdditional = [];
experimentalSettings(1).data_modalities = {'surfaces','thicknesses', 'volumes'};
experimentalSettings(1).connectionsFlag = false;
experimentalSettings(1).dichotomizeFlag = false;
experimentalSettings(1).indx_conn = [];
experimentalSettings(1).nperm = 5000;
experimentalSettings(1).nboot = 5000;
experimentalSettings(1).caxisLimList = [0.04 0.04 0.04];
experimentalSettings(1).nullmodel = 'subjectLabelPermutation';
experimentalSettings(1).phenotype = 'default';

experimentalSettings(end+1).experimentName = 'main_analyses';
experimentalSettings(end).timePointImaging = '2'; % Imaging related data
experimentalSettings(end).timePoints = [1 3]; % Averaging over behavior related data
experimentalSettings(end).covariatesAdditional = [];
experimentalSettings(end).data_modalities = {'connvS'};
experimentalSettings(end).connectionsFlag = true;
experimentalSettings(end).dichotomizeFlag = false;
experimentalSettings(end).nperm = 5000;
experimentalSettings(end).indx_conn = indx_conn_SC;
experimentalSettings(end).nboot = 5000;
experimentalSettings(end).caxisLimList = [0.02];
experimentalSettings(end).nullmodel = 'subjectLabelPermutation';

experimentalSettings(end+1).experimentName = 'main_analyses';
experimentalSettings(end).timePointImaging = '2'; % Imaging related data
experimentalSettings(end).timePoints = [1 3]; % Averaging over behavior related data
experimentalSettings(end).covariatesAdditional = [];
experimentalSettings(end).data_modalities = {'connvF'};
experimentalSettings(end).indx_conn = indx_conn_FC;
experimentalSettings(end).connectionsFlag = true;
experimentalSettings(end).nperm = 5000;
experimentalSettings(end).dichotomizeFlag = false;
experimentalSettings(end).nboot = 5000;
experimentalSettings(end).caxisLimList = [0.02];
experimentalSettings(end).nullmodel = 'subjectLabelPermutation';

experimentalSettings(end+1).experimentName = 'main_analyses';
experimentalSettings(end).timePointImaging = '2'; % Imaging related data
experimentalSettings(end).timePoints = [1 3]; % Averaging over behavior related data
experimentalSettings(end).covariatesAdditional = [];
experimentalSettings(end).data_modalities = {'tfmri'};
experimentalSettings(end).connectionsFlag = false;
experimentalSettings(end).nperm = 5000;
experimentalSettings(end).dichotomizeFlag = false;
experimentalSettings(end).nboot = 5000;

%% Use spin-based permutation testing
experimentalSettings(end+1) = experimentalSettings(1);
experimentalSettings(end).nullmodel = 'spintesting';

experimentalSettings(end+1) = experimentalSettings(2);
experimentalSettings(end).nullmodel = 'spintesting';

experimentalSettings(end+1) = experimentalSettings(3);
experimentalSettings(end).nullmodel = 'spintesting';
 
%% Use data from second-imaging visit
experimentalSettings(end+1).experimentName = 'repeated_scans';
experimentalSettings(end).timePointImaging = '3'; % Imaging related data
experimentalSettings(end).timePoints = [4]; % Averaging over behavior related data
experimentalSettings(end).covariatesAdditional = [];
experimentalSettings(end).data_modalities = {'thicknessesT3', 'surfacesT3', 'volumesT3'};
experimentalSettings(end).connectionsFlag = false;
experimentalSettings(end).dichotomizeFlag = false;
experimentalSettings(end).nperm = 5000;
experimentalSettings(end).nboot = 5000;
experimentalSettings(end).caxisLimList = [0.02 0.04 0.04];

experimentalSettings(end+1).experimentName = 'repeated_scans';
experimentalSettings(end).timePointImaging = '3'; % Imaging related data
experimentalSettings(end).timePoints = [4]; % Averaging over behavior related data
experimentalSettings(end).covariatesAdditional = [];
experimentalSettings(end).data_modalities = {'connvST3'};
experimentalSettings(end).connectionsFlag = true;
experimentalSettings(end).dichotomizeFlag = false;
experimentalSettings(end).nperm = 5000;
experimentalSettings(end).indx_conn = indx_conn_SC;
experimentalSettings(end).nboot = 5000;
experimentalSettings(end).caxisLimList = [0.02];

experimentalSettings(end+1).experimentName = 'repeated_scans';
experimentalSettings(end).timePointImaging = '3'; % Imaging related data
experimentalSettings(end).timePoints = [4]; % Averaging over behavior related data
experimentalSettings(end).covariatesAdditional = [];
experimentalSettings(end).data_modalities = {'connvFT3'};
experimentalSettings(end).indx_conn = indx_conn_FC;
experimentalSettings(end).connectionsFlag = true;
experimentalSettings(end).nperm = 5000;
experimentalSettings(end).dichotomizeFlag = false;
experimentalSettings(end).nboot = 5000;
experimentalSettings(end).caxisLimList = [0.02];

experimentalSettings(end+1).experimentName = 'repeated_scans';
experimentalSettings(end).timePointImaging = '3'; % Imaging related data
experimentalSettings(end).timePoints = [4]; % Averaging over behavior related data
experimentalSettings(end).covariatesAdditional = [];
experimentalSettings(end).data_modalities = {'tfmriT3'};
experimentalSettings(end).connectionsFlag = false;
experimentalSettings(end).nperm = 5000;
experimentalSettings(end).dichotomizeFlag = false;
experimentalSettings(end).nboot = 5000;

%% Use participants from the hold-out dataset
% Important note: To compute the data on the hold-out dataset, one has to
% update the `results_version` variable to end with "_holdout" and reload
% all data again using script_main_initialization.m

experimentalSettings(1).experimentName = 'holdout';
experimentalSettings(end).timePointImaging = '2'; % Imaging related data
experimentalSettings(end).timePoints = [1 3]; % Averaging over behavior related data
experimentalSettings(end).covariatesAdditional = [];
experimentalSettings(end).data_modalities = {'thicknesses', 'surfaces', 'volumes'};
experimentalSettings(end).connectionsFlag = false;
experimentalSettings(end).dichotomizeFlag = false;
experimentalSettings(end).nperm = 5000;
experimentalSettings(end).nboot = 5000;
experimentalSettings(end).caxisLimList = [0.02 0.04 0.04];

experimentalSettings(end+1).experimentName = 'holdout';
experimentalSettings(end).timePointImaging = '2'; % Imaging related data
experimentalSettings(end).timePoints = [1 3]; % Averaging over behavior related data
experimentalSettings(end).covariatesAdditional = [];
experimentalSettings(end).data_modalities = {'connvS'};
experimentalSettings(end).connectionsFlag = true;
experimentalSettings(end).dichotomizeFlag = false;
experimentalSettings(end).nperm = 5000;
experimentalSettings(end).indx_conn = indx_conn_SC;
experimentalSettings(end).nboot = 5000;
experimentalSettings(end).caxisLimList = [0.02];

experimentalSettings(end+1).experimentName = 'holdout';
experimentalSettings(end).timePointImaging = '2'; % Imaging related data
experimentalSettings(end).timePoints = [1 3]; % Averaging over behavior related data
experimentalSettings(end).covariatesAdditional = [];
experimentalSettings(end).data_modalities = {'connvF'};
experimentalSettings(end).indx_conn = indx_conn_FC;
experimentalSettings(end).connectionsFlag = true;
experimentalSettings(end).dichotomizeFlag = false;
experimentalSettings(end).nperm = 5000;
experimentalSettings(end).nboot = 5000;
experimentalSettings(end).caxisLimList = [0.02];

experimentalSettings(end+1).experimentName = 'holdout';
experimentalSettings(end).timePointImaging = '2'; % Imaging related data
experimentalSettings(end).timePoints = [1 3]; % Averaging over behavior related data
experimentalSettings(end).covariatesAdditional = [];
experimentalSettings(end).data_modalities = {'tfmri'};
experimentalSettings(end).connectionsFlag = false;
experimentalSettings(end).nperm = 5000;
experimentalSettings(end).dichotomizeFlag = false;
experimentalSettings(end).nboot = 5000;

%% Alternative symptom severity score defintion

experimentalSettings(end+1).experimentName = 'alt-pheno';
experimentalSettings(end).timePointImaging = '2'; % Imaging related data
experimentalSettings(end).timePoints = [1 3]; % Averaging over behavior related data
experimentalSettings(end).covariatesAdditional = [];
experimentalSettings(end).data_modalities = {'surfaces','thicknesses', 'volumes'};
experimentalSettings(end).connectionsFlag = false;
experimentalSettings(end).dichotomizeFlag = false;
experimentalSettings(end).indx_conn = [];
experimentalSettings(end).nperm = 5000;
experimentalSettings(end).nboot = 1;
experimentalSettings(end).caxisLimList = [0.04 0.04 0.04];
experimentalSettings(end).nullmodel = 'subjectLabelPermutation';
experimentalSettings(end).phenotype = 'alt-pheno';

experimentalSettings(end+1).experimentName = 'alt-pheno';
experimentalSettings(end).timePointImaging = '2'; % Imaging related data
experimentalSettings(end).timePoints = [1 3]; % Averaging over behavior related data
experimentalSettings(end).covariatesAdditional = [];
experimentalSettings(end).data_modalities = {'connvS'};
experimentalSettings(end).connectionsFlag = true;
experimentalSettings(end).dichotomizeFlag = false;
experimentalSettings(end).nperm = 5000;
experimentalSettings(end).indx_conn = indx_conn_SC;
experimentalSettings(end).nboot = 1;
experimentalSettings(end).caxisLimList = [0.02];
experimentalSettings(end).nullmodel = 'subjectLabelPermutation';
experimentalSettings(end).phenotype = 'alt-pheno';

experimentalSettings(end+1).experimentName = 'alt-pheno';
experimentalSettings(end).timePointImaging = '2'; % Imaging related data
experimentalSettings(end).timePoints = [1 3]; % Averaging over behavior related data
experimentalSettings(end).covariatesAdditional = [];
experimentalSettings(end).data_modalities = {'connvF'};
experimentalSettings(end).indx_conn = indx_conn_FC;
experimentalSettings(end).connectionsFlag = true;
experimentalSettings(end).nperm = 5000;
experimentalSettings(end).dichotomizeFlag = false;
experimentalSettings(end).nboot = 1;
experimentalSettings(end).caxisLimList = [0.02];
experimentalSettings(end).nullmodel = 'subjectLabelPermutation';
experimentalSettings(end).phenotype = 'alt-pheno';

%% Exclude relatives
% Important note: script_exclude_relatives changes the brain data, so you
% need to reload the data if you want to do other analyses on the whole
% dataset.

experimentalSettings(end+1) = experimentalSettings(1);
experimentalSettings(end).nboot = 1;
experimentalSettings(end).excludeRelatives = true;

experimentalSettings(end+1) = experimentalSettings(2);
experimentalSettings(end).nboot = 1;
experimentalSettings(end).excludeRelatives = true;

experimentalSettings(end+1) = experimentalSettings(3);
experimentalSettings(end).nboot = 1;
experimentalSettings(end).excludeRelatives = true;

experimentalSettings(end+1) = experimentalSettings(4);
experimentalSettings(end).nboot = 1;
experimentalSettings(end).excludeRelatives = true;

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

    if excludeRelatives
        script_exclude_relatives
    end

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
    
    % Run analyses
    script_calculate_regional_maps_IDA

    close all
    
    % Present results
    if connectionsFlag

        data_modality = data_modalities{1};
        caxis_lim = caxisLimList;
        script_present_associations_connectivity

    elseif any(contains(data_modalities, 'tfmri'))

        script_present_associations_tfmri

    else

        script_present_associations_structural

    end
    
    catch ME
        warning(jsonencode(ME));
    end
    
end
