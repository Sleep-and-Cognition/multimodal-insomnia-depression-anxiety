delete(fullfile(results_path_general, 'load_prepare_data.txt'));
diary(fullfile(results_path_general, 'load_prepare_data.txt'));

%% Load Imaging data
fprintf('load data...');

rall = load(region_properties_file);
conn = load(sc_connectivity_file);
connFC = load(fc_connectivity_file);

% Find intersection of subjects present in rall, conn and connFC
subjectID = intersect(rall.info_table.subjectID, conn.info_table.subjectID);
subjectID = intersect(subjectID, connFC.info_table.subjectID);

% Prepare regional data
[~, ~, J] = intersect(subjectID, rall.info_table.subjectID, 'stable');
rall.info_table = rall.info_table(J,:);
rall.regionProperties = rall.regionProperties(:, :, J);

% Prepare structural connectivity 
[~, ~, J] = intersect(subjectID, conn.info_table.subjectID, 'stable');
conn.info_table = conn.info_table(J,:);
conn.connectivity = conn.connectivity(:, J);

% Prepare functional connectivity 
[~, ~, J] = intersect(subjectID, connFC.info_table.subjectID, 'stable');
connFC.info_table = connFC.info_table(J,:);
connFC.connectivity = connFC.connectivity(:, J);

% Prepare tfMRI data
tfmri_table = load(tfmri_file);
tfmri_table = tfmri_table.info_table;
[~, ~, J] = intersect(subjectID, tfmri_table.subjectID, 'stable');
assert(isequal(tfmri_table.subjectID(J), subjectID));
tfmri = tfmri_table.f_25052_2_0(J)'; % Median BOLD effect (in group-defined amygdala activation mask)

% Construct derived variables
volumes = squeeze(rall.regionProperties(:, 4, :));
surfaces = squeeze(rall.regionProperties(:, 5, :));
thicknesses = squeeze(rall.regionProperties(:, 6, :));

regionProperties = rall.regionProperties;
nregions = length(rall.regionDescriptions);
connv = conn.connectivity;
connv_FC = connFC.connectivity;

% Add covariates table
covar_table = load(covariates_file);
covar_table = covar_table.info_table;

[~, I, J] = intersect(subjectID, covar_table.Properties.RowNames, 'stable');
assert(isequal(covar_table.Properties.RowNames(J), subjectID));
covar_table = covar_table(J,:);

fprintf(' succeeded \n');

%% Filter subjects
fprintf('Filter data...\n');

% 26500: Whether the T2-FLAIR was used or not (in addition to the T1) to
% run FreeSurfer
indx_include_1 = covar_table.f_26500_2_0 == 'Yes';

fprintf('T2-FLAIR was used in addition to T1 to run FS\n');
fprintf('\tSubjects included: %i. Subjects excluded %i.\n', ...
    nnz(indx_include_1), nnz(~indx_include_1));


% 25926: Intensity scaling for T2_FLAIR
indx_include_2 = covar_table.f_25926_2_0 > 3;

fprintf('Intensity scaling for T2 Flair > 3\n');
fprintf('\tSubjects included: %i. Subjects excluded %i.\n', ...
    nnz(indx_include_2), nnz(~indx_include_2));


% NaNs in FC or SC matrix or volumes, surfaces, thicknesses
indx_include3 = all(~isnan(connv_FC),1) ...
                & all(~isnan(connv),1) ...
                & all(~isnan(surfaces(15:end,:)),1) ...
                & all(~isnan(thicknesses(15:end,:)),1) ...
                & all(~isnan(volumes(1:14,:)),1);

indx_include3 = indx_include3';

% Missing 25758 Scanner longitudinal (Z) brain position
indx_include_4 = ~isnan(covar_table.f_25758_2_0);

fprintf('Mising Scanner longitudinal (Z) brain position\n');
fprintf('\tSubjects included: %i. Subjects excluded %i.\n', ...
    nnz(indx_include_4), nnz(~indx_include_4));

fprintf(['Exclude subjects with NaNs in FC or ' ...
    'SC matrix or volumes, surfaces, thicknesses.\n']);
fprintf('\tSubjects included: %i. Subjects excluded %i.\n', ...
    nnz(indx_include3), nnz(~indx_include3));

indx_include = indx_include_1 & indx_include_2 ...
    & indx_include3 & indx_include_4;

fprintf('Total: Subjects included: %i. Subjects excluded %i.\n', ...
    nnz(indx_include), nnz(~indx_include));

connv(:,~indx_include) = NaN;
connv_FC(:,~indx_include) = NaN;
volumes(:,~indx_include) = NaN;
surfaces(:,~indx_include) = NaN;
thicknesses(:,~indx_include) = NaN;
tfmri(~indx_include) = NaN;

fprintf('Filter data... succeeded \n');

%% Quality-control
fprintf('Quality control...\n');

outliers_combined = quality_control( ...
    'results_path_setup', results_path_setup, ...
'phenotype_qc_file', phenotype_qc_file, ...
'subjectID', subjectID, ...
'figures_path_setup', figures_path_setup, ...
'connv_FC', connv_FC, ...
'connv', connv, ...
'covar_table', covar_table, ...
'phenotype_qc_file', phenotype_qc_file, ...
'volumes', volumes, ...
'surfaces', surfaces, ...
'thicknesses', thicknesses, ...
'tfmri', tfmri, ...
'indx_include', indx_include);
close all

assert(nnz(outliers_combined & ~indx_include) == 0);

connv(:, outliers_combined) = NaN;
connv_FC(:, outliers_combined) = NaN;
volumes(:, outliers_combined) = NaN;
surfaces(:, outliers_combined) = NaN;
thicknesses(:, outliers_combined) = NaN;
tfmri(outliers_combined) = NaN;

fprintf('Quality control... succeeded \n');

%% Phenotype
fprintf('Make phenotype...\n');

script_make_phenotype

fprintf('Make phenotype... succeeded \n');

diary('off');