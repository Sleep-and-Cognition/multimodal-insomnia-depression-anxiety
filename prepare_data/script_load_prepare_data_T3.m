%% Load imaging data from repeated visit (T3)

fprintf('Load all data...');

%% Regional measures
rallT3 = load(region_properties_T3_file);
[~, I, J] = intersect(subjectID, rallT3.info_table.subjectID, 'stable');
assert(isequal(subjectID(I), rallT3.info_table.subjectID(J)))

thicknessesT3 = nan(size(thicknesses));
thicknessesT3(:, I) = squeeze(rallT3.regionProperties(:, 6,J));

indx = contains(rallT3.propertyDescriptions, 'surface');
surfacesT3 = nan(size(surfaces));
surfacesT3(:, I) = squeeze(rallT3.regionProperties(:, indx,J));

indx = contains(rallT3.propertyDescriptions, 'volume');
volumesT3 = nan(size(volumes));
volumesT3(:, I) = squeeze(rallT3.regionProperties(:, indx,J));

%% SC and FC repeated
fcT3 = load(fc_connectivity_T3_file);

[~, I, J] = intersect(subjectID, fcT3.info_table.subjectID, 'stable');
connv_FC_T3 = nan(size(connv_FC,1), length(subjectID));
connv_FC_T3(:, I) = fcT3.connectivity(:,J);

scT3 = load(sc_connectivity_T3_file);
[~, I, J] = intersect(subjectID, scT3.info_table.subjectID, 'stable');
connv_T3 = nan(size(connv,1), length(subjectID));
connv_T3(:, I) = squeeze(scT3.connectivity(:, 3,J));

clear fcT3 scT3


%% TFMRI data
tfmri_table = load(tfmri_file);
tfmri_table = tfmri_table.info_table;
[~, ~, J] = intersect(subjectID, tfmri_table.subjectID, 'stable');
assert(isequal(tfmri_table.subjectID(J), subjectID));
tfmriT3 = tfmri_table.f_25052_3_0(J)'; %Median BOLD effect (in group-defined amygdala activation mask)

clear tfmri_table
fprintf(' succeeded.\n');

%% Filter subjects
fprintf('Filter data...\n');

% 26500: Whether the T2-FLAIR was used or not (in addition to the T1) to
% run FreeSurfer
indx_include_1 = covar_table.f_26500_3_0 == 'Yes';

fprintf('T2-FLAIR was used in addition to T1 to run FS\n');
fprintf('\tSubjects included: %i. Subjects excluded %i.\n', ...
    nnz(indx_include_1), nnz(~indx_include_1));


% 25926: Intensity scaling for T2_FLAIR
indx_include_2 = covar_table.f_25926_3_0 > 3;

fprintf('Intensity scaling for T2 Flair > 3\n');
fprintf('\tSubjects included: %i. Subjects excluded %i.\n', ...
    nnz(indx_include_2), nnz(~indx_include_2));

% NaNs in FC or SC matrix or volumes, surfaces, thicknesses
indx_include3 = all(~isnan(connv_FC_T3),1) ...
                & all(~isnan(connv_T3),1) ...
                & all(~isnan(surfacesT3(15:end,:)),1) ...
                & all(~isnan(thicknessesT3(15:end,:)),1) ...
                & all(~isnan(volumesT3(1:14,:)),1);

indx_include3 = indx_include3';

% Missing 25758 Scanner longitudinal (Z) brain position
indx_include_4 = ~isnan(covar_table.f_25758_3_0);

fprintf('Mising Scanner longitudinal (Z) brain position\n');
fprintf('\tSubjects included: %i. Subjects excluded %i.\n', ...
    nnz(indx_include_4), nnz(~indx_include_4));

% Exclude 25926
% In the main analysis is the Intensity scaling of T2 FLAIR (ID: 25926) a
% covariate. However, in the repeated sample only 1 subject has this item
% deviating from the main value. Making it better to exclude this subject.
indx_include_5 = covar_table.f_25926_3_0 == 6;

fprintf('Intensity scaling of T2 FLAIR == 6\n');
fprintf('\tSubjects included: %i. Subjects excluded %i.\n', ...
    nnz(indx_include_5), nnz(~indx_include_5));

% Similarly only 19 people have different values for the table position is
% the Z-coordinate (ID: 25759).
indx_include_6 = covar_table.f_25759_3_0 == -1042;
fprintf('The table position is the Z-coordinate == -1042\n');
fprintf('\tSubjects included: %i. Subjects excluded %i.\n', ...
    nnz(indx_include_6), nnz(~indx_include_6));


fprintf(['Exclude subjects with NaNs in FC ' ...
    'or SC matrix or volumes, surfaces, thicknesses.\n']);
fprintf('\tSubjects included: %i. Subjects excluded %i.\n', ...
    nnz(indx_include3), nnz(~indx_include3));

indx_include = indx_include_1 & indx_include_2 ...
    & indx_include3 & indx_include_4 &indx_include_5 & indx_include_6;

fprintf('Total: Subjects included: %i. Subjects excluded %i.\n', ...
    nnz(indx_include), nnz(~indx_include));

connv_T3(:,~indx_include) = NaN;
connv_FC_T3(:,~indx_include) = NaN;
volumesT3(:,~indx_include) = NaN;
surfacesT3(:,~indx_include) = NaN;
thicknessesT3(:,~indx_include) = NaN;
tfmriT3(~indx_include) = NaN;

fprintf('Filter data... succeeded \n');


%% Quality-control

% save results and figures in specific repeated measures data.
checkdir(fullfile(results_path_setup, 'T3'));
checkdir(fullfile(figures_path_setup, 'T3'));

outliers_combined = quality_control( ...
    'results_path_setup', fullfile(results_path_setup, 'T3'), ...
'phenotype_qc_file', phenotype_qc_file, ...
'subjectID', subjectID, ...
'figures_path_setup', fullfile(figures_path_setup, 'T3'), ...
'connv_FC', connv_FC_T3, ...
'connv', connv_T3, ...
'covar_table', covar_table, ...
'phenotype_qc_file', phenotype_qc_file, ...
'volumes', volumesT3, ...
'surfaces', surfacesT3, ...
'thicknesses', thicknessesT3, ...
'tfmri', tfmriT3, ...
'indx_include', indx_include);

close all

assert(nnz(outliers_combined & ~indx_include) == 0);

connv_T3(:, outliers_combined) = NaN;
connv_FC_T3(:, outliers_combined) = NaN;
volumesT3(:, outliers_combined) = NaN;
surfacesT3(:, outliers_combined) = NaN;
thicknessesT3(:, outliers_combined) = NaN;
tfmriT3(outliers_combined) = NaN;

fprintf('Quality control... succeeded \n');

%% Create derived variables
fprintf('Make derived variables (connvFT3 and connvST3)...');

connvFT3 = connv_FC_T3(indx_conn_FC, :);

connvST3 = connv_T3(indx_conn_SC, :);
connvST3(connvST3 == 0) = NaN;

fprintf(' succeeded.\n');

diary('off');