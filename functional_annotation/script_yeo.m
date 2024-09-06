%% Compute overlap between brain maps and functional networks

%% Initialize
delete(fullfile(results_path, 'yeo.txt'))
diary(fullfile(results_path, 'yeo.txt'));

yeo = load(yeo_file);

[~, I,J] = intersect(regionDescriptions, yeo.regionDescriptions);
A = nan(size(regionDescriptions,1), size(yeo.matrixDK2Yeo,2));
A(I,:) = yeo.matrixDK2Yeo(J,:);
yeo.matrixDK2Yeo = A;

if all(isnan(A))
    fprintf('No yeo data available. Abort mission.\n');
    return
end

% Remove NONE module
yeo.matrixDK2Yeo = yeo.matrixDK2Yeo(:, 2:end);
yeo.yeoDescriptions =yeo.yeoDescriptions(2:end); 
nModules = size(yeo.matrixDK2Yeo,2);


b_difference_t_bt = difference_t_bt;
b_difference_t_bt = normalize(b_difference_t_bt, 'zscore');

%% Prepare feature maps
clear feature_map_all

for im = 1:nModules
        if connectionsFlag
            feature_map = yeo.matrixDK2Yeo(:, im);
            feature_map = (feature_map > prctile(feature_map, 85));

            feature_map = squareform(setdiag(feature_map .* feature_map', 0))';
            feature_map(~indx_conn) = 0;
            feature_map_all(:, im) = feature_map ./ nansum(feature_map);

        else
            feature_map = yeo.matrixDK2Yeo(:, im);
            feature_map = (feature_map > prctile(feature_map, 85));
            feature_map_all(:, im) = feature_map ./ nansum(feature_map);
        end
end

%% Calculate prev_pr and permutations

% Initialize variables
cr = nan(nModules, 1);
cr_lo = nan(nModules, 1);
cr_up = nan(nModules, 1);
cp = nan(nModules, 1);

crall = nan(nModules, nsymptoms);
cpall = nan(nModules, nsymptoms);
cpall_uncorrected = nan(nModules, nsymptoms);

cr_bt = nan(nModules, nboot);
crbtall = nan(nModules, nboot, nsymptoms);

for isymptom = 1:nsymptoms
            
    for im = 1:nModules
        feature_map = feature_map_all(:,im);

        cr(im) = nansum(difference_t(:,isymptom) .* feature_map);
        cr_bt(im, :) = nansum(squeeze(b_difference_t_bt(:, isymptom,:)) .* feature_map);
        
        % Get p-values from permutation testing        
        r_perm = nansum(squeeze(difference_t_perm(:, isymptom,:)) .* feature_map);
        [~,cp(im) ] = ttest2(r_perm, cr(im));
        
    end
    
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
        num2cell(cp(J)) yeo.yeoDescriptions(J)']);
    
    fprintf('Number significant %i\n', nnz(cp < 0.05));

    J = J(cp(J) < 0.05);
    if length(J) >0
        J = flipud(J);
        J = J(1:min(3, length(J)));

        res = [yeo.yeoDescriptions(J)' num2cell(cr(J)) num2cell(cp(J))]';
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
pres = nan(3,3,nModules);

for im = 1:nModules

    A = crbtall(im, :, 1) - crbtall(im, :, 2);
    pres(1,2,im) = 2*normcdf(-abs(mean(A) ./ std(A)));
    
    
    A = crbtall(im, :, 3) - crbtall(im, :, 2);
    pres(2,3,im) = 2*normcdf(-abs(mean(A) ./ std(A)));
    
    
    A = crbtall(im, :, 1) - crbtall(im, :, 3);
    pres(1,3,im) = 2*normcdf(-abs(mean(A) ./ std(A)));
    
end

[~, ~, pres_corr] = fdr(pres(:));
pres_corr = reshape(pres_corr, size(pres));

% Report results
fprintf('Corrected across %i tests.\n', nnz(~isnan(pres(:))));
fprintf('Significant differences between features:\n');

for im = 1:nModules
    fprintf('%s:\n', yeo.yeoDescriptions{im});
    for is = 1:nsymptoms
        for iis = 1:nsymptoms
            if pres_corr(is, iis, im) < 0.05

                fprintf('\tSIGNIFICANT: %s (%.3f) - %s (%.3f) p = %.3f\n', ...
                    symptoms{is}, mean(crbtall(im, :, is),2), ...
                    symptoms{iis}, mean(crbtall(im, :, iis),2), ...
                    pres_corr(is,iis, im));
            end
        end
    end
end

yeoDescriptions = yeo.yeoDescriptions;

save(fullfile(results_path, 'results_yeo'), ...
    'crall', 'cpall', 'cpall_uncorrected', ...
    'yeoDescriptions', 'symptoms', 'pres_corr');

%% Present results in a figure

% Sort modules by the module categories
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
for im = 1:numel(crall)
  
    y = u(im) + [-0.5 -0.5 0.5 0.5];
    A = v(im) + [-0.5 0.5 0.5 -0.5];
    h(im) = patch(A, y, crall(im), ...
        'EdgeColor', [1 1 1], ...
        'LineWidth', 0.5);
end

% Add hatches to non-significant items
for im = 1:length(h)
    if cpall(im) >= 0.05
        hatch(h(im))
    end
end


for i = 1:nc
   
    for j = 1:nsymptoms
        t = '';
        for jj = 1:nsymptoms
            if (cpall(i,j) < 0.05) && ((pres_corr(j, jj, i) < 0.05) || ...
                    (pres_corr(jj, j, i) < 0.05)) ...
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
yticklabels(yeo.yeoDescriptions);

% Set title
title('Association between brain map with functional networks')

% Add colorbar
c = colorbar('SouthOutside');
c.Box = false;
c.TickLength = 0;
c.Label.String = 'Correlation (r)';

g = gcf;
g.Units = 'normalized';
g.Position = [    0.5344    0.5795    0.1350    0.2992];

% Save figure to file
pause(1)
print(fullfile(figures_path, ['yeo_' resType]), '-dsvg');

diary('off');