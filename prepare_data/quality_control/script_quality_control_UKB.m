%% Quality control script for UKB measures
% Outlier Detection and Visualization of Quality Control Measures

% Input Variables:
% - phenotype_qc_file: A file containing the quality control (QC) measures for phenotypes
% - subjectID: A vector containing subject IDs
% - covar_table: A table containing covariate information for each subject
% - IQR_RATIO_THRESHOLD: A numeric threshold for determining outliers based on interquartile range (IQR)
% - figures_path_setup: Location where output figures are saved

% Output Variables:
% - excluded_UKB: A logical array indicating whether each subject is considered an outlier or not

% Explanation:
% This MATLAB script determines outliers based on the interquartile range
% (IQR) of UKB provided QC measures. The output variables include a logical
% array indicating whether each subject is considered an outlier and the
% total number of excluded subjects.
%
% Notes:
% The brain x,y,z position UKB items are included because they are
% covariates and one of them shows some clear outliers. The table position
% variable is not included because it had a small peak of non-outliers
% which includes resulted in 900 (false-)outliers.

%% Load data
fprintf('Quality control script with UKB provided items.\n');

% Assert all variables are set
assert(exist('phenotype_qc_file','var') == 1);
assert(isfile(phenotype_qc_file));
assert(exist('subjectID','var') == 1);

qc_table = load(phenotype_qc_file);
qc_table = qc_table.info_table;

% Match QC table with subject list.
[~, I, J] = intersect(subjectID, qc_table.Properties.RowNames, 'stable');
assert(isequal(qc_table.Properties.RowNames(J), subjectID));
qc_table = qc_table(J,:);

% Combine QC table (excluding PC items) with covar-table (that includes the
% brain position fields).
qc_table = [qc_table(:,10:end) covar_table];

indx = ismember(qc_table.Properties.VariableNames, ...
    {'f_25731_2_0', 'f_25732_2_0', 'f_25733_2_0', 'f_25739_2_0', ...
    'f_25737_2_0', 'f_25744_2_0', 'f_25741_2_0', 'f_25746_2_0', 'f_25756_2_0', ...
    'f_25757_2_0', 'f_25758_2_0'});

assert(nnz(indx) == 11, 'Not all QC variables found.');

qc_measures = qc_table{:, indx};
qc_names = qc_table.Properties.VariableDescriptions(indx);

fprintf('The following QC measures were included:\n');
disp(qc_names')

%% Determine outliers

[S, l_thres, u_thres] = isoutlier(qc_measures, 'quartiles', ...
    'ThresholdFactor', IQR_RATIO_THRESHOLD);
excluded_UKB = any(S,2);

%% Visualize outliers
figure('color', 'white');

sn = ceil(sqrt(size(qc_measures,2)));
sm = ceil(size(qc_measures,2) ./ sn);

for i = 1:size(qc_measures,2)
    subplot(sn, sm, i);
    
    h = histogram(qc_measures(:,i));

    hold on;
    h2 = histogram(qc_measures(S(:,i),i), 'BinEdges', h.BinEdges);
       
    xlabel(textwrap(qc_names(i), 75));
    
    title(sprintf('N = %i', sum(S(:,i))));

    % Simplify figure
    a = gca;
    a.Box = 'off';
    
    % Make it visually more appealing
    a.XAxis.LineWidth = 1;
    a.YAxis.LineWidth = 1;
    a.FontSize = 12;
    
    % Give the bars more importance
    color1 = [72 150 236] ./ 255;
    h.LineWidth = 1;
    h.EdgeColor = color1;
    h.FaceAlpha = 0;
    h.FaceColor = color1;
    
    color2 = [0.3765    0.2510    0.6902];
    h2.LineWidth = 1;
    h2.EdgeColor = color2;
    h2.FaceAlpha = 0.5;
    h2.FaceColor = color2;
    
    line([l_thres(i) l_thres(i)], ylim, 'Color', color2, 'linewidth', 1);
    line([u_thres(i) u_thres(i)], ylim, 'Color', color2, 'linewidth', 1);

    
    
end

f = gcf;
f.Units = 'normalized';
f.Position = [0 0 0.95 0.95];

print(fullfile(figures_path_setup, 'QC', 'UKB_measures'), '-dsvg');


%% Visualize number of violations per subject
s = sum(S, 2);

figure('color', 'white');
h = histogram(s(s>0));
xlabel('number of violations per subject (larger than 0)');
a = gca;
a.Box = 'off';

% Make it visually more appealing
a.XAxis.LineWidth = 1;
a.YAxis.LineWidth = 1;
a.FontSize = 12;

% Give the bars more importance
color1 = [72 150 236] ./ 255;
h.LineWidth = 1;
h.EdgeColor = color1;
h.FaceAlpha = 0.3;
h.FaceColor = color1;

print(fullfile(figures_path_setup, 'QC', 'UKB_nViolations'), '-dsvg');

fprintf('percentage outliers (based on UKB phenotypes: %.3g%%\n', ...
    100*nnz(any(S, 2)) ./ size(S ,1));


