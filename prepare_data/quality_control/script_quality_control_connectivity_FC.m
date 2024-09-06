%% Quality control script for FC connectivity matrices
% Outlier Detection of Connectivity Measures
%
% Input Variables:
% - connv_FC: A matrix containing connectivity values for each subject
% - IQR_RATIO_THRESHOLD: A numeric threshold for determining outliers based on interquartile range (IQR)
% - figures_path_setup: Figures will be saved in the subdirectory `QC` 
%
% Output Variables:
% - excluded_connv: A logical array indicating whether each subject is considered an outlier or not
% - figures: Two figures visualizing the connectivity measures and the number of violations per subject
%
% Explanation: This MATLAB script calculates outlier measures for each
% subject based on connectivity values. It then determines outliers based
% on the interquartile range (IQR) and a specified threshold factor. The
% script visualizes the connectivity measures as histograms and the number
% of violations per subject. The output variables include a logical array
% indicating whether each subject is considered an outlier and the total
% number of excluded subjects.

%% initialize
fprintf('Quality control script with FC connectivity.\n');

nwei = 1;
nsub = size(connv_FC,2); % number of subjects

%% Calculate outlier measures

C = zeros(nsub, nwei); % table with the score for each subject
metricDescription = cell(0);

C(:,1) = nzmean(connv_FC, 1);
metricDescription{1} = 'weighting';

C(all(isnan(connv_FC),1), :) = nan;

fprintf('The following measures were included:\n');
disp(metricDescription');

%% determine outliers

[S, l_thres, u_thres] = isoutlier(C, 'quartiles', ...
    'ThresholdFactor', IQR_RATIO_THRESHOLD);
excluded_connv_FC = any(S,2);
n_excluded = nnz(excluded_connv_FC);

%% Visualize

figure('color', 'white');
for i = 1:size(C,2)
    subplot(1, size(C,2), i)
    h = histogram(C(:,i));

    hold on;
    h2 = histogram(C(S(:,i),i), 'BinEdges', h.BinEdges);
       
    xlabel(metricDescription{i});
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

print(fullfile(figures_path_setup, 'QC', 'connv_FC_measures'), '-dsvg');


%% Visualize number of violations per subject

Ss = sum(S,2);
figure('color', 'white');
h = histogram(Ss(Ss>0));
xlabel('number of violations per subject (larger than 0)');

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
h.FaceAlpha = 0.3;
h.FaceColor = color1;

print(fullfile(figures_path_setup, 'QC', 'connv_FC_nViolations'), '-dsvg');

fprintf('percentage outliers (based on Connectivity: %.3g%%\n', ...
    100*nnz(excluded_connv_FC) ./ length(excluded_connv_FC));