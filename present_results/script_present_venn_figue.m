%% Create visualization of regions associated with symptom-specific and shared findings
%% Initialization

% PREPARE: Run script_main_initialization.m
regionDescriptions = hcp.regionDescriptions;
regionDescriptions_all = regionDescriptions;

% PREPARE: run first sections of NKE and receptor script
for ir = 1 : size(receptor_data,2)
feature_map = receptor_data(:,ir);
feature_map = (feature_map > prctile(feature_map, 75));
receptor_all(:, ir)  = feature_map ./ nansum(feature_map);
end

for indx = 1:size(mean_feature_per_ROI, 2)
feature_map = mean_feature_per_ROI(:,indx);
nke_all(:,indx) = feature_map ./ nansum(feature_map);
end

receptor_descriptions  = string(receptor_descriptions);

% Colors
cm_venn.I = [140	203	210	] ./ 255;
cm_venn.D = [235	160	146	]./ 255;
cm_venn.A = [241	186	131	]./ 255;
cm_venn.ID = [188	138	127	]./ 255;
cm_venn.IA = [134	175	160	]./ 255;
cm_venn.DA = [229	132	96	]./ 255;
cm_venn.IDA = [187	128	108	]./ 255;

%% Create figures
% Each figure illustrates the localization of symptom-specific and shared
% findings across modalities. The findings are also presented in Figure 4a 

%% Insomnia symptom specific findings: subcortical reward
indx = find(strcmp(feature_description, 'Reward'));
f_map = mean_feature_per_ROI(:, indx);
f_map(15:end) = 0;
C = f_map;

cm = interpolate_cbrewer([1 1 1; cm_venn.I]*255, 'spline', 1000) ./ 255;
C(1:7) = mean([C(1:7), C(8:14)],2);
C(8:14) = C(1:7);
C(15:48) = mean([C(15:48), C(49:82)],2);
C(49:82) = C(15:48);
plotBrain(regionDescriptions, C, cm, 'atlas', 'aparc_aseg');

%% Depressive symptom specific findings: cortical language and reward
indx = find(strcmp(feature_description, 'Language'));
f_map = mean_feature_per_ROI(:, indx);
f_map(1:14) = 0;
C = f_map;
indx = find(strcmp(feature_description, 'Reward'));
f_map = mean_feature_per_ROI(:, indx);
f_map(1:14) = 0;
C = max(C, f_map);

cm = interpolate_cbrewer([1 1 1; cm_venn.D]*255, 'spline', 1000) ./ 255;
C(1:7) = mean([C(1:7), C(8:14)],2);
C(8:14) = C(1:7);
C(15:48) = mean([C(15:48), C(49:82)],2);
C(49:82) = C(15:48);
plotBrain(regionDescriptions, C, cm, 'atlas', 'aparc_aseg');

%% Anxiety symptom specific findings: D1, DAT, mGlueR5, H3, amygdala
indx = find(strcmp(receptor_descriptions, 'D1'));
f_map = receptor_all(:, indx);
C = f_map;

indx = find(strcmp(receptor_descriptions, 'DAT'));
f_map = receptor_all(:, indx);
C = C + f_map;

indx = find(strcmp(receptor_descriptions, 'mGluR5'));
f_map = receptor_all(:, indx);
C = C + f_map;

indx = find(strcmp(receptor_descriptions, 'H3'));
f_map = receptor_all(:, indx);
C = C + f_map;

C = C ./ max(C(:));
indx = find(contains(regionDescriptions, 'Amygdala'));
C(indx) = 1;

C(isnan(C)) = 0;
cm = interpolate_cbrewer([1 1 1; cm_venn.A]*255, 'spline', 1000) ./ 255;
C(1:7) = mean([C(1:7), C(8:14)],2);
C(8:14) = C(1:7);
C(15:48) = mean([C(15:48), C(49:82)],2);
C(49:82) = C(15:48);
plotBrain(regionDescriptions, C, cm, 'atlas', 'aparc_aseg');

%% Ins + Anx specific findings: H3, A4B2
indx = find(strcmp(receptor_descriptions, 'H3'));
f_map = receptor_all(:, indx);
C = f_map;

indx = find(strcmp(receptor_descriptions, 'A4B2'));
f_map = receptor_all(:, indx);
C = C + f_map;

C(isnan(C)) = 0;
cm = interpolate_cbrewer([1 1 1; cm_venn.IA]*255, 'spline', 1000) ./ 255;
C(1:7) = mean([C(1:7), C(8:14)],2);
C(8:14) = C(1:7);
C(15:48) = mean([C(15:48), C(49:82)],2);
C(49:82) = C(15:48);
plotBrain(regionDescriptions, C, cm, 'atlas', 'aparc_aseg');

%% Dep + Anx specific findings: D2
indx = find(strcmp(receptor_descriptions, 'D2'));
f_map = receptor_all(:, indx);
C = f_map;

C(isnan(C)) = 0;
cm = interpolate_cbrewer([1 1 1; cm_venn.DA]*255, 'spline', 1000) ./ 255;
C(1:7) = mean([C(1:7), C(8:14)],2);
C(8:14) = C(1:7);
C(15:48) = mean([C(15:48), C(49:82)],2);
C(49:82) = C(15:48);
plotBrain(regionDescriptions, C, cm, 'atlas', 'aparc_aseg');

%% Regional findings common across Ins + Dep + Anx symptoms: thalamus
C = zeros(size(regionDescriptions));
indx = find(contains(regionDescriptions, 'Thalamus-Proper'));
C(indx) = 7;

cm = interpolate_cbrewer([1 1 1; cm_venn.IDA]*255, 'spline', 1000) ./ 255;
C(1:7) = mean([C(1:7), C(8:14)],2);
C(8:14) = C(1:7);
C(15:48) = mean([C(15:48), C(49:82)],2);
C(49:82) = C(15:48);
plotBrain(regionDescriptions, C, cm, 'atlas', 'aparc_aseg');







