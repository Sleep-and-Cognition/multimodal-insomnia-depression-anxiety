function outliers_combined = quality_control(varargin)
%% Quality control
%
% Input:
% results_path_setup    Location where results are written
% phenotype_qc_file     MAT file with QC phenotypes
% subjectID             Subject list with subject IDs from UKB
% figures_path_setup    Location where output figures are saved
% connv_FC				A matrix containing connectivity values for each subject
% connv					A matrix containing connectivity values for each subject
% covar_table			A table containing covariate information for each subject
% phenotype_qc_file		A file containing the quality control (QC) measures for phenotypes
% volumes				A matrix containing volume measures for each subject and region
% surfaces				A matrix containing surface area measures for each subject and region
% thicknesses			A matrix containing cortical thickness measures for each subject and region
% tfmri                 tfmri data
%
% Output:
% outliers_combined

% Parse and validate optional arguments
while ~isempty(varargin)
    if numel(varargin) == 1
        error('quality_control:missingOptions', ...
            'Optional arguments must come in pairs.');
    end
    
    switch lower(varargin{1})
        case 'results_path_setup'
            results_path_setup = varargin{2};
        case 'phenotype_qc_file'
            phenotype_qc_file = varargin{2};
        case 'subjectid'
            subjectID = varargin{2};
        case 'figures_path_setup'
            figures_path_setup = varargin{2};
        case 'connv_fc'
            connv_FC = varargin{2};
        case 'connv'
            connv = varargin{2};
        case 'covar_table'
            covar_table = varargin{2};
        case 'volumes'
            volumes = varargin{2};
        case 'surfaces'
            surfaces = varargin{2};
        case 'thicknesses'
            thicknesses = varargin{2};
        case 'tfmri'
            tfmri = varargin{2};
        case 'indx_include'
            indx_include = varargin{2};
        otherwise
            error('quality_control:unknownOption', ...
                'The option %s is unknown.', varargin{1});
    end
    
    % Remove processed arguments
    varargin(1:2) = [];
end

diary(fullfile(results_path_setup, 'QC.txt'));
checkdir(fullfile(figures_path_setup, 'QC'));

IQR_RATIO_THRESHOLD = 2.5;
fprintf('IQR ratio applied: %g\n', IQR_RATIO_THRESHOLD);

checkdir(fullfile(figures_path_setup, 'QC'));

script_quality_control_UKB % provides excluded_UKB 
script_quality_control_structural % provides excluded_regions 
script_quality_control_connectivity_SC % provides excluded_connv 
script_quality_control_connectivity_FC % provides excluded_connv_FC 

excluded_tfmri = isoutlier(tfmri, 'quartiles', ...
    'ThresholdFactor', IQR_RATIO_THRESHOLD)';

% Subjects that are already excluded in the filter step should not be
% included here.
excluded_UKB(~indx_include) = 0;
excluded_regions(~indx_include) = 0;
excluded_connv(~indx_include) = 0;
excluded_connv_FC(~indx_include) = 0;
excluded_tfmri(~indx_include) = 0;

outliers_combined = excluded_UKB | excluded_regions | ...
excluded_connv | excluded_connv_FC | excluded_tfmri;
fprintf('percentage outliers (combined): %.3g%%\n', ...
    100*nnz(outliers_combined) ./ size(outliers_combined ,1));

fprintf('#outliers UKB: %i\n', nnz(excluded_UKB));
fprintf('#outliers regions: %i\n', nnz(excluded_regions));
fprintf('#outliers connections SC: %i\n', nnz(excluded_connv));
fprintf('#outliers connections FC: %i\n', nnz(excluded_connv_FC));
fprintf('#outliers tfmri: %i\n', nnz(excluded_tfmri));
fprintf('#outliers combined: %i\n', nnz(outliers_combined));


diary off