%% Quality control of brain regions to detect and exclude outliers
% This MATLAB script determines outliers based on the interquartile range
% (IQR) of regional data (volumes, surfaces, and thicknesses) and a
% specified threshold factor.
%
% Input Variables:
% - volumes: A matrix containing volume measures for each subject and region
% - surfaces: A matrix containing surface area measures for each subject and region
% - thicknesses: A matrix containing cortical thickness measures for each subject and region
% - IQR_RATIO_THRESHOLD: A numeric threshold for determining outliers based on interquartile range (IQR)
% - figures_path_setup: A path to save the output figures

% Output Variables:
% - excluded_regions: A logical array indicating whether each subject is considered an outlier or not

%% Calculate outlier measures
fprintf('Quality control script with regional measures (volumes, surfaces and thicknesses).\n');

% Based on the inter quartile range (IQR).

S = false(size(volumes, 2), 3);
S(:,1) = any(isoutlier(volumes', 'quartiles', 'ThresholdFactor', IQR_RATIO_THRESHOLD), 2);
S(:,2) = any(isoutlier(surfaces', 'quartiles', 'ThresholdFactor', IQR_RATIO_THRESHOLD), 2);
S(:,3) = any(isoutlier(thicknesses', 'quartiles', 'ThresholdFactor', IQR_RATIO_THRESHOLD), 2);

excluded_regions = any(S,2);

%% Visualize
cm = cbrewer('seq', 'Blues', 1000);

figure('color', 'white');
data_modalities = {'volumes', 'surfaces', 'thicknesses'};
for im = 1:length(data_modalities)
    
    eval(['this_data = ' data_modalities{im} ';']);
    
    hc = nan(100, size(this_data, 1));
    for i = 1:size(this_data, 1)
        data_scale =  linspace(min(this_data(i,:)), max(this_data(i,:)), 101);
        if all(isnan(data_scale))
            hc(:,i) = nan;
        else
        hc(:, i) = histcounts(this_data(i,:), data_scale);
        end
    end
    
    subplot(length(data_modalities), 1, im);
    imagesc(1:size(this_data, 1), data_scale, hc)
    colormap(cm);
    title(data_modalities{im} );
    
end

print(fullfile(figures_path_setup, 'QC', 'regions_heatmap'), '-dsvg');

%% Visualize number of violations per subject

s = sum(S,2);
figure('color', 'white');
h = histogram(s(s>0));
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

print(fullfile(figures_path_setup, 'QC', 'regions_nViolations'), '-dsvg');

fprintf('percentage outliers (based on regions: %.3g%%\n', ...
    100*nnz(excluded_regions) ./ length(excluded_regions));
