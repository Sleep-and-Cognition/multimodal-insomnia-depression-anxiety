%% Make phenotype

%% Region descriptions
regionDescriptions_all = rall.regionDescriptions;

%% Phenotype data
diary(fullfile(results_path_setup, 'make_phenotype.txt'));
checkdir(fullfile(figures_path_setup, 'make_phenotype'));

tmp = load(phenotypes_file);
info_table = tmp.info_table;
clear tmp

%% Select only the subjects in the imaging dataset

[~, I, J] = intersect(subjectID, info_table.subjectID, 'stable');
assert(isequal(info_table.subjectID(J), subjectID));
info_table = info_table(J,:);

%% Make group mask for FC and SC.

indx_conn_FC = nanmean(connv_FC,2)>0.3;
disp(nnz(indx_conn_FC) ./ numel(indx_conn_FC))
connvF = connv_FC(indx_conn_FC, :);

indx_conn_SC = mean(connv>0,2) > 0.6;
disp(nnz(indx_conn_SC) ./ numel(indx_conn_SC))
connvS = connv(indx_conn_SC, :);
connvS(connvS == 0) = NaN;

%% Prepare phenotype data
n12_names = {'1930', 'Miserableness'; ...
'1940', 'Irritability'; ...
'1950', 'Sensitivity / hurt feelings'; ...
'1960', 'Fed-up feelings'; ...
'1970', 'Nervous feelings'; ...
'1980', 'Worrier / anxious feelings'; ...
'1990', 'Tense / highly strung'; ...
'2000', 'Worry too long after embarrassment'; ...
'2010', 'Suffer from nerves'; ...
'2020', 'Loneliness, isolation'; ...
'2030', 'Guilty feelings'; ...
'2040', 'Risk taking'};

n12 = nan(size(info_table,1), size(n12_names,1),4);
for it = 1:4
    indx = find(contains(info_table.Properties.VariableNames, ...
        strcat('f_', n12_names(:,1), '_', num2str(it-1), '_0')));
    for iq = 1:nnz(indx)
        tmp = info_table{:, indx(iq)};
        n12(tmp == 'Yes', iq, it) = 1;
        n12(tmp == 'No', iq, it) = 0;
    end
end

rds4_names = {'2050', 'Frequency of depressed mood in last 2 weeks'; ...
'2060', 'Frequency of unenthusiasm / disinterest in last 2 weeks'; ...
'2070', 'Frequency of tenseness / restlessness in last 2 weeks'; ...
'2080', 'Frequency of tiredness / lethargy in last 2 weeks'};

rds4 = nan(size(info_table,1), size(rds4_names,1),4);

for it = 1:4
    indx = find(contains(info_table.Properties.VariableNames, ...
        strcat('f_', rds4_names(:,1), '_', num2str(it-1), '_0')));
    for iq = 1:nnz(indx)
        tmp = info_table{:, indx(iq)};
        rds4(tmp == 'Not at all', iq, it) = 0;
        rds4(tmp == 'Several days', iq, it) = 1;
        rds4(tmp == 'More than half the days', iq, it) = 2;
        rds4(tmp == 'Nearly every day', iq, it) = 3;
    end
end

% Insomnia
insomnia_items = [info_table.f_1200_0_0 info_table.f_1200_1_0 
    info_table.f_1200_2_0 info_table.f_1200_3_0];
insomnia_missing = (insomnia_items == 'Prefer not to answer') | ...
ismissing(insomnia_items);
insomnia_score = nan(size(insomnia_items));
insomnia_score(insomnia_items == 'Usually') = 2;
insomnia_score(insomnia_items == 'Sometimes') = 1;
insomnia_score(insomnia_items == 'Never/rarely') = 0;
insomnia_score(insomnia_missing) = missing;

% Derived measures
n12_nervous = squeeze(sum(n12(:, [5 7 9], :),2));
phq2 = squeeze(sum(rds4(:, [1 2], :), 2));

%% Add information about medication.

addpath('phenotype/');
load('phenotypes_medicin.mat');

medicin_table.subjectID = medicin_table.Properties.RowNames;
[~, I, J] = intersect(subjectID, medicin_table.subjectID, 'stable');
assert(isequal(medicin_table.subjectID(J), subjectID));

medicin_table = medicin_table(J,:);

mc = readtable('medicine_classes.csv');
mc.Group = categorical(mc.Group);

x = medicin_table{:,1:end-1};
clear res
mc_a = mc(mc.Group == 'Anxiolytics',:);
res(:, 1) = any(contains(string(x), mc_a.Drug_Name, 'IgnoreCase', true), 2);

mc_a = mc(mc.Group == 'Antidepressants',:);
res(:, 2) = any(contains(string(x), mc_a.Drug_Name, 'IgnoreCase', true), 2);

mc_a = mc(mc.Group == 'Antipsychotics',:);
res(:, 3) = any(contains(string(x), mc_a.Drug_Name, 'IgnoreCase', true), 2);

mc_a = mc(mc.Group == 'Mood_Stabiliser',:);
res(:, 4) = any(contains(string(x), mc_a.Drug_Name, 'IgnoreCase', true), 2);

medication = res;
medication_names = {'Anxiolytics', 'Antidepressants', 
    'Antipsychotics', 'Mood_Stabiliser'};

%% Social economic status
x = readtable('phenotypes_covariates-20240712-170231.csv_output.csv');
info_cov_table = x;

names = info_cov_table.Properties.VariableNames;
names = strrep(names, 'p', 'f_');
names = strrep(names, '_i', '_');
names = strcat(names, '_0_0');

names = strrep(names, '_0_0_0', '_0_0');
names = strrep(names, '_1_0_0', '_1_0');
names = strrep(names, '_2_0_0', '_2_0');
names = strrep(names, '_3_0_0', '_3_0');

info_cov_table.Properties.VariableNames = names;

info_cov_table.subjectID = cellfun(@num2str, ...
    num2cell(info_cov_table.eid_0_0), 'uniformOutput', false);
[~, I, J] = intersect(subjectID, info_cov_table.subjectID, 'stable');
assert(isequal(info_cov_table.subjectID(J), subjectID(I)));

townsend = nan(size(subjectID));
townsend(I) = info_cov_table.f_22189_0_0(J);

%% Childhood trauma
% From: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9170771/#appsec1
%
% Never True (0)	
% Rarely True (1)	
% Sometimes True (2)	
% Often True (3)	
% Very Often True (4)
%
% Field ID	Question/Statement	Reverse Coding Used?
% 20487	Felt hated by a family member as a child	No	
% 20488	Physically abused by family as a child	No	
% 20489	Felt loved as a child	Yes	
% 20490	Sexually molested as a child	No
% 20491	Someone to take to doctor when needed as a child	Yes	

CTS_items = nan(size(info_cov_table,1),5);

% First non-reversed items:
items = categorical([info_cov_table.f_20487_0_0 ...
    info_cov_table.f_20488_0_0 info_cov_table.f_20490_0_0]);

for iq = 1:size(items,2)
    tmp = items(:, iq);
    CTS_items(tmp == 'Never true', iq) = 0;
    CTS_items(tmp == 'Often', iq) = 3;
    CTS_items(tmp == 'Rarely true', iq) = 1;
    CTS_items(tmp == 'Sometimes true', iq) = 2;
    CTS_items(tmp == 'Very often true', iq) = 4;
end

% Reversed items:
items = categorical([info_cov_table.f_20489_0_0 info_cov_table.f_20491_0_0]);

for iq = 1:size(items,2)
    tmp = items(:, iq);
    CTS_items(tmp == 'Never true', iq+3) = 4;
    CTS_items(tmp == 'Often', iq+3) = 1;
    CTS_items(tmp == 'Rarely true', iq+3) = 3;
    CTS_items(tmp == 'Sometimes true', iq+3) = 2;
    CTS_items(tmp == 'Very often true', iq+3) = 0;
end

CTS_items = sum(CTS_items,2);
CTS = nan(size(subjectID));
CTS(I) = CTS_items(J);


%% Adult trauma
% From: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9170771/#appsec1
%
% Never True (0)	
% Rarely True (1)	
% Sometimes True (2)	
% Often True (3)	
% Very Often True (4)
%
% Field ID	Question/Statement	Reverse Coding Used?
% 20521	Belittlement by partner or ex-partner as an adult	No	
% 20522	Been in a confiding relationship as an adult	Yes	
% 20523	Physical violence by a partner or ex-partner as an adult	No	
% 20524	Sexual interference by a partner or ex-partner without consent as an adult	No
% 20525	Able to pay rent/mortgage as an adult	Yes	

ATS_items = nan(size(info_cov_table,1),5);

% First non-reversed items:
items = categorical([info_cov_table.f_20521_0_0 info_cov_table.f_20523_0_0 info_cov_table.f_20524_0_0]);

for iq = 1:size(items,2)
    tmp = items(:, iq);
    ATS_items(tmp == 'Never true', iq) = 0;
    ATS_items(tmp == 'Often', iq) = 3;
    ATS_items(tmp == 'Rarely true', iq) = 1;
    ATS_items(tmp == 'Sometimes true', iq) = 2;
    ATS_items(tmp == 'Very often true', iq) = 4;
end

% Reversed items:
items = categorical([info_cov_table.f_20525_0_0 info_cov_table.f_20522_0_0]);

for iq = 1:size(items,2)
    tmp = items(:, iq);
    ATS_items(tmp == 'Never true', iq+3) = 4;
    ATS_items(tmp == 'Often', iq+3) = 1;
    ATS_items(tmp == 'Rarely true', iq+3) = 3;
    ATS_items(tmp == 'Sometimes true', iq+3) = 2;
    ATS_items(tmp == 'Very often true', iq+3) = 0;
end

ATS_items = sum(ATS_items,2);
ATS = nan(size(subjectID));
ATS(I) = ATS_items(J);

TS = ATS + CTS;

%% Smoking
smoking_items = info_cov_table.f_20162_0_0;
info_cov_table.f_20116_0_0 = categorical(info_cov_table.f_20116_0_0);

smoking_items(:) = 1;
smoking_items(info_cov_table.f_20116_0_0 == 'Never') = 0;

smoking = nan(size(subjectID));
smoking(I) = smoking_items(J);

%% Diabetes
% 130706 Date E10 first reported (insulin-dependent diabetes mellitus)
% 130707 Source of report of E10 (insulin-dependent diabetes mellitus)
% 130708 Date E11 first reported (non-insulin-dependent diabetes mellitus)
% 130709 Source of report of E11 (non-insulin-dependent diabetes mellitus)
% 130710 Date E12 first reported (malnutrition-related diabetes mellitus)
% 130711 Source of report of E12 (malnutrition-related diabetes mellitus)
% 130712 Date E13 first reported (other specified diabetes mellitus)
% 130713 Source of report of E13 (other specified diabetes mellitus)
% 130714 Date E14 first reported (unspecified diabetes mellitus)
% 130715 Source of report of E14 (unspecified diabetes mellitus)
diabetes_items = [categorical(info_cov_table.f_130707_0_0), ...
categorical(info_cov_table.f_130709_0_0), ...
categorical(info_cov_table.f_130711_0_0), ...
categorical(info_cov_table.f_130713_0_0), ...
categorical(info_cov_table.f_130715_0_0)];
diabetes_items = any(~ismissing(diabetes_items),2);

diabetes = nan(size(subjectID));
diabetes(I) = diabetes_items(J);

%% Hypertension
% 131286 Date I10 first reported (essential (primary) hypertension)
% 131287 Source of report of I10 (essential (primary) hypertension)
% 131288 Date I11 first reported (hypertensive heart disease)
% 131289 Source of report of I11 (hypertensive heart disease)
% 131290 Date I12 first reported (hypertensive renal disease)
% 131291 Source of report of I12 (hypertensive renal disease)
% 131292 Date I13 first reported (hypertensive heart and renal disease)
% 131293 Source of report of I13 (hypertensive heart and renal disease)
% 131294 Date I15 first reported (secondary hypertension)
% 131295 Source of report of I15 (secondary hypertension)

hypertension_items = [categorical(info_cov_table.f_131287_0_0), ...
    categorical(info_cov_table.f_131289_0_0), ...
    categorical(info_cov_table.f_131291_0_0), ...
    categorical(info_cov_table.f_131293_0_0), ...
    categorical(info_cov_table.f_131295_0_0)];

hypertension_items = any(~ismissing(hypertension_items),2);

hypertension = nan(size(subjectID));
hypertension(I) = hypertension_items(J);

%% Maximum education
education_items = nan(size(info_cov_table,1),1);
tmp = (info_cov_table.f_6138_0_0);

% Order of the following lines is important because the search terms are
% subtrings of eachother.
education_items(contains(tmp, 'None of the above')) = 1;
education_items(contains(tmp, 'CSEs or equivalent')) = 2;
education_items(contains(tmp, 'O levels/GCSEs or equivalent')) = 3; 
education_items(contains(tmp, 'NVQ or HND or HNC or equivalent')) = 3.5;
education_items(contains(tmp, 'A levels/AS levels or equivalent')) = 4;
education_items(contains(tmp, 'College or University degree')) = 5;
education_items(contains(tmp, 'Other professional qualifications eg: nursing, teaching')) = 6;

education = nan(size(subjectID));
education(I) = education_items(J);

% Years of education: Data-Field 845
educationy = info_cov_table.f_845_0_0;
educationy = educationy - 5; % UK compulsory school start at age 5
educationy(info_cov_table.f_845_0_0 <= 0) = NaN;

% impute missing years of education based on qualification level
% using mean of non-missing data from same qualification level.

educationy(education_items == 1) = NaN;
educationy(ismissing(educationy) & education_items == 2) = 11.1;
educationy(ismissing(educationy) & education_items == 3) = 11.5;
educationy(ismissing(educationy) & education_items == 3) = 12;
educationy(ismissing(educationy) & education_items == 4) = 13.1;
educationy(ismissing(educationy) & education_items == 5) = 17;
educationy(ismissing(educationy) & education_items == 6) = 12.9;

education_years = nan(size(subjectID));
education_years(I) = educationy(J);

%% Alcohol
% A log-transformed measure of grams of ethanol typically consumed per day,
% either the max across time or the individual time points (0.0, 1.0, 2.0).
x = readtable(fullfile(absolute_path, 'data', ...
    'alc_measures_siemon_20240712.txt'));
info_cov_table = x;

names = info_cov_table.Properties.VariableNames;
names = strrep(names, 'p', 'f_');
names = strrep(names, '_i', '_');
names = strcat(names, '_0_0');

names = strrep(names, '_0_0_0', '_0_0');
names = strrep(names, '_1_0_0', '_1_0');
names = strrep(names, '_2_0_0', '_2_0');
names = strrep(names, '_3_0_0', '_3_0');

info_cov_table.Properties.VariableNames = names;

info_cov_table.subjectID = cellfun(@num2str, ...
    num2cell(info_cov_table.f_eid_0_0), 'uniformOutput', false);
[~, I, J] = intersect(subjectID, info_cov_table.subjectID, 'stable');
assert(isequal(info_cov_table.subjectID(J), subjectID(I)));

alcohol = nan(size(subjectID));
alcohol(I) = info_cov_table.gm_eth_ln_0_0(J);