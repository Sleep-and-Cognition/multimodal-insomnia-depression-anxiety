%% Compute overlap between brain maps and neurotransmission systems

%% Initialization

delete(fullfile(results_path, 'receptor.txt'))
diary(fullfile(results_path, 'receptor.txt'));

rd = load(fullfile(absolute_path, ...
    ['nke_data/results/receptor_data/receptorData_' atlas '.mat']));
receptor_data = rd.receptorData;

[~, I] = ismember(regionDescriptions, regionDescriptions_all);
receptor_data = receptor_data(I,:);

if all(isnan(receptor_data))
    fprintf('No receptor data available. Abort mission.\n');
    return
end

receptor_descriptions = rd.receptorDescriptions;
receptor_categories = rd.receptorCategories;

[~, J] = sort(receptor_categories);
receptor_descriptions = receptor_descriptions(J);
receptor_categories = receptor_categories(J);
receptor_data = receptor_data(:,J);

nreceptors = size(receptor_data, 2);

b_difference_t_bt = difference_t_bt;
b_difference_t_bt = normalize(b_difference_t_bt, 'zscore');


%% Normalize and prepare feature maps

clear feature_map_all

for ir = 1:nreceptors
        if connectionsFlag
            feature_map = receptor_data(:,ir);
            feature_map = (feature_map > prctile(feature_map, 75));
            feature_map = squareform(setdiag(feature_map .* feature_map', 0))';
            feature_map(~indx_conn) = 0;
            feature_map_all(:, ir) = feature_map ./ nansum(feature_map);

        else
            feature_map = receptor_data(:,ir);
            feature_map = (feature_map > prctile(feature_map, 75));
            feature_map_all(:, ir) = feature_map ./ nansum(feature_map);
        end
end

%% Calculate prev_pr and permutations

% Initialize variables
cr = nan(nreceptors, 1);
cr_lo = nan(nreceptors, 1);
cr_up = nan(nreceptors, 1);
cp = nan(nreceptors, 1);

crall = nan(nreceptors, nsymptoms);
cpall = nan(nreceptors, nsymptoms);
cpall_uncorrected = nan(nreceptors, nsymptoms);

cr_bt = nan(nreceptors, nboot);
crbtall = nan(nreceptors, nboot, nsymptoms);

for isymptom = 1:nsymptoms
            
    % Correlation between brain map and #topics involved per component
    for ir = 1:nreceptors
        feature_map = feature_map_all(:,ir);
        
        cr(ir) = nansum(difference_t(:,isymptom) .* feature_map);
        cr_bt(ir, :) = nansum(squeeze(b_difference_t_bt(:, isymptom,:)) .* ...
            feature_map);
        
        % Get p-values from permutation testing        
        r_perm = nansum(squeeze(difference_t_perm(:, isymptom,:)) .* ...
            feature_map);
        [~, cp(ir)] = ttest2(r_perm, cr(ir));

    end
    
    % The p-values are corrected for the number of receptors using FDR.
    cpall_uncorrected(:, isymptom) = cp;
    crall(:, isymptom) = cr;
    crbtall(:, :, isymptom) = cr_bt;

end

% Correct p-values for multiple comparisons.
indx = crall ~= 0;
fprintf('Corrected across %i tests.\n', nnz(indx));
[~, ~, cpall(indx)] = fdr(cpall_uncorrected(indx));

%% Report results
for isymptom = 1:nsymptoms

    fprintf('\n--------------\nsymptom %s:\n', symptoms{isymptom});
    cr = crall(:,isymptom);
    cp = cpall(:,isymptom);    

    [~, J] = sort(cr, 'descend', 'MissingPlacement', 'last');
       
    fprintf('Correlation between prev_pr and comp_region\n');
    disp([{'r'} {'r05'} {'r95'} {'pFDR'} {'comp'}; ...
        num2cell(cr(J)) num2cell(cr_lo(J))  num2cell(cr_up(J)) ...
        num2cell(cp(J)) receptor_descriptions(J)']);
    
    fprintf('Number significant %i\n', nnz(cp < 0.05));

    J = J(cp(J) < 0.05);
    if length(J) > 0
        J = flipud(J);
        J = J(1:min(3, length(J)));

        res = [receptor_descriptions(J)' num2cell(cr(J)) num2cell(cp(J))]';
        res = sprintf('"%s" (β = %.3f, p = %.3f), ', res{:});
        res = strrep(res, 'p = 0.000', 'p < 0.001');
        res = lower(res);
        disp(res);
    end

end



%% Perform comparative analyses between symptom scores
% Compare the observed correlations between the different symptoms (i.e.
% symptom scores). This is done by comparing whether the difference of
% their estimated distribution (obtained using bootstrapping) is
% significantly different from zero.

pres = nan(3,3,nreceptors);

for ir = 1:nreceptors

    x = crbtall(ir, :, 1) - crbtall(ir, :, 2);
    pres(1,2,ir) = 2*normcdf(-abs(mean(x) ./ std(x)));
    
    
    x = crbtall(ir, :, 3) - crbtall(ir, :, 2);
    pres(2,3,ir) = 2*normcdf(-abs(mean(x) ./ std(x)));
    
    
    x = crbtall(ir, :, 1) - crbtall(ir, :, 3);
    pres(1,3,ir) = 2*normcdf(-abs(mean(x) ./ std(x)));
    
end

[~, ~, pres_corr] = fdr(pres(:));
pres_corr = reshape(pres_corr, size(pres));

% Report results
fprintf('Corrected across %i tests.\n', nnz(~isnan(pres(:))));
fprintf('Significant differences between features:\n');

for ir = 1:nreceptors
    fprintf('%s:\n', receptor_descriptions{ir});
    for is = 1:nsymptoms
        for iis = 1:nsymptoms
            if pres_corr(is, iis, ir) < 0.05
                fprintf('\tSIGNIFICANT: %s (%.3f) - %s (%.3f) p = %.3f\n', ...
                    symptoms{is}, mean(crbtall(ir, :, is),2), ...
                    symptoms{iis}, mean(crbtall(ir, :, iis),2), ...
                    pres_corr(is,iis, ir));
            end
        end
    end
end

%% Save all statistics
save(fullfile(results_path, 'results_receptor'), ...
    'crall', 'cpall', 'cpall_uncorrected', ...
    'receptor_descriptions', 'symptoms', 'pres_corr');

%% Present results in a figure

% Sort receptors by the receptor categories
nc = size(crall, 1);

rall = squeeze(mean(crbtall,2));
rall = rall ./ max(abs(rall));
rall = max(0, -rall);

% Create a heatmap of the correlation array
figure('color', 'white');
colormap(cm);
hold on;

% x and y location of each item in the heatmap
u = repmat((1:nc)', 1, nsymptoms);
v = repmat(1:nsymptoms, nc, 1);

% Plot the heatmap
clear h
for ir = 1:numel(crall)
  
    y = u(ir) + [-0.5 -0.5 0.5 0.5];
    x = v(ir) + [-0.5 0.5 0.5 -0.5];
    h(ir) = patch(x, y, crall(ir), ...
        'EdgeColor', [1 1 1], ...
        'LineWidth', 0.5);
end

% Add hatches to non-significant items
for ir = 1:length(h)
    if cpall(ir) >= 0.05
        hatch(h(ir))
    end
end


for i = 1:nc
   
    for j = 1:nsymptoms
        t = '';
        for jj = 1:nsymptoms
            if (cpall(i,j) < 0.05) && ((pres_corr(j, jj, i) < 0.05) ...
                    || (pres_corr(jj, j, i) < 0.05)) ...
                    && (crall(i, j) < crall(i, jj))
                t = [t '▾'];
            end
        end
        text(j, i, t, 'FontSize', 12, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'color', 'white');
    end
end


% Set axes
caxis([-caxis_lim caxis_lim]);
axis tight

xticks(1:nsymptoms);
xticklabels(symptoms);
xtickangle(45);
yticks(1:nc);
yticklabels(receptor_descriptions);

% Set title
title('Correlation symptom map with receptor density')

% Add colorbar
c = colorbar('SouthOutside');
c.Box = false;
c.TickLength = 0;
c.Label.String = 'Correlation (r)';

g = gcf;
g.Units = 'normalized';
g.Position = [    0.5344    0.6479    0.0746    0.2806];

% Save figure to file
pause(1)
print(fullfile(figures_path, ['receptor_grid_004' resType]), '-dsvg');

diary('off');