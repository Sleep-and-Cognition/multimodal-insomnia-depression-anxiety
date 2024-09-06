%% Initialize
fprintf('Initialize...');

results_version = '240606';
atlas = 'aparc';

load('initialization.mat');

if contains(results_version, 'holdout')
    region_properties_file = holdout.region_properties_file;
    sc_connectivity_file = holdout.sc_connectivity_file;
    fc_connectivity_file = holdout.fc_connectivity_file;
end

region_properties_file = strrep(region_properties_file, 'ATLAS', atlas);
neurosynth_file = strrep(neurosynth_file, 'ATLAS', atlas);
yeo_file = strrep(yeo_file, 'ATLAS', atlas);

sc_connectivity_file = strrep(sc_connectivity_file, 'ATLAS', atlas);
fc_connectivity_file = strrep(fc_connectivity_file, 'ATLAS', atlas);

scripts_path = fullfile(absolute_path, 'scripts');
phenotypes_file = fullfile(scripts_path, 'phenotype', ...
    'phenotypes_REMOVE_augmented.mat');
covariates_file = fullfile(scripts_path, 'preanalysis', ...
    'covariates', 'covariates.mat');

cd(scripts_path);
addpath(scripts_path);
addpath(genpath(scripts_path));

results_path_general = fullfile(absolute_path, 'results', ...
    results_version, atlas);
checkdir(results_path_general);

figures_path_general = fullfile(absolute_path, 'figures', ...
    results_version, atlas);
checkdir(figures_path_general);

results_path_setup = results_path_general; 
checkdir(results_path_setup)
figures_path_setup = figures_path_general; 
checkdir(figures_path_setup)

assert(isequal(which('fdr'), fullfile(scripts_path, 'misc', 'fdr.m')));

cm = flipud(cbrewer('div', 'RdBu', 1000));
cm_blue = cbrewer('seq', 'Blues', 1000);

cmmetal = [66  83 175; ...
72 150 236; ...
74 168 238; ...
84 186 209; ...
65 148 136; ...
103 172  91; ...
151 192  92; ...
208 218  89; ...
252 234  96; ...
246 194  68; ...
241 156  56; ...
236  99  55; ...
116  86  74; ...
158 158 158; ...
102 124  13; ...
225  82  65; ...
214  57 100; ...
143  54 170; ...
96  64 176] ./ 255;

fprintf(' succeeded \n');

%% Run tests
% Test the getResults implementation:
runtests

%% Load brain atlas information

data = load('regionData.mat');
regionDescriptions = data.regionDescriptions;
regionDescriptionsAll = regionDescriptions;
regionCoordinates = data.coordinates;

%% Load all data 

script_load_prepare_data
script_load_prepare_data_T3

%% Analyses
% The scripts below are sometimes mutually exclusive and can best be run
% separately.

%% Compute and present brain-behavior associations

script_main_run_experiments
script_main_run_experiments_holdout

%% Report other stuff

% Create figure with localization of effects for all the symptom specific
% and symptom shared annotations.
script_present_venn_figue

% Calculate the correlation of the association maps across modalities.
script_present_correlation_across_modalities

% Estimate the clinical impact of the observed association strengths.
script_present_implications

% Create overview table with all brain-behavior associations.
script_present_overview_table

% Validate the symptom severity scores with respect to in-depth measures
% and diagnostic data
script_validate_symptom_scores

%% Sensitivity analyses

% Validate the findings from the main analysis in the hold-out imaging
% dataset, second imaging visit dataset and with an alternative phenotype
% using bootstrapped results from the main analysis.
script_analyse_bootstrap

% Examine the effect of excluding relatives
script_present_results_excrelatives

% Examine the effect of adding additional covariates in the analyses
script_present_correlation_covariates
script_analyse_extra_covariates








