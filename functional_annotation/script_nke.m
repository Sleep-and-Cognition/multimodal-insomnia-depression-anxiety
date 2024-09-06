%% Compute overlap between brain maps and cognitive emotional domains

%% Initialization

delete(fullfile(results_path, 'nke.txt'));
diary(fullfile(results_path, 'nke.txt'));

nfeatures = 8;
nke = load(fullfile(absolute_path, ...
    ['nke_data/results/result_make_nke_', atlas, ...
    '_211117_', num2str(nfeatures), '.mat']));

info = nke.info;
mean_feature_per_ROI = nke.mean_feature_per_ROI;

[~, I] = ismember(regionDescriptions, nke.regionDescriptions);
mean_feature_per_ROI = mean_feature_per_ROI(I,:);
[~, I] = intersect(info.CLUSTER, 1:nfeatures);
feature_description = info.TITLE(I);

[~, J] = sort(feature_description);
mean_feature_per_ROI = mean_feature_per_ROI(:,J);
feature_description = feature_description(J);

b_difference_t_bt = difference_t_bt;
b_difference_t_bt = normalize(b_difference_t_bt, 'zscore');

%% Normalize and prepare feature maps
clear feature_map_all

for i = 1:size(mean_feature_per_ROI, 2)

    if connectionsFlag
        feature_map = mean_feature_per_ROI(:, i);
        feature_map = squareform(setdiag(feature_map * feature_map', 0))';
        feature_map(~indx_conn) = 0;
        feature_map_all(:,i) = feature_map ./ nansum(feature_map);
    else
        feature_map = mean_feature_per_ROI(:,i);
        feature_map_all(:,i) = feature_map ./ nansum(feature_map);
    end

end

%% Calculate prev_pr and permutations

% Initialize variables
cr = nan(nfeatures, 1);
cr_lo = nan(nfeatures, 1);
cr_up = nan(nfeatures, 1);
cp = nan(nfeatures, 1);

crall = nan(nfeatures, nsymptoms);
cpall = nan(nfeatures, nsymptoms);
cpall_uncorrected = nan(nfeatures, nsymptoms);

cr_bt = nan(nfeatures, nboot);
crbtall = nan(nfeatures, nboot, nsymptoms);

for isymptom = 1:nsymptoms
   
    nregions = size(mean_feature_per_ROI,1);
    
    for i = 1:size(mean_feature_per_ROI, 2)
                
        feature_map = feature_map_all(:,i);

        cr(i) = nansum(difference_t(:,isymptom) .* feature_map);
        cr_bt(i, :) = nansum(squeeze(b_difference_t_bt(:, isymptom,:)) .* ...
            feature_map);
        
        % Get p-values from permutation testing        
        r_perm = nansum(squeeze(difference_t_perm(:, isymptom,:)) .* ...
            feature_map);
        [~,cp(i) ] = ttest2(r_perm, cr(i));

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
    fprintf('\n--------------\nsymptom %s:\n', symptoms{isymptom})
    cr = crall(:,isymptom);
    cp = cpall(:,isymptom);
    
    [~, J] = sort(cr, 'descend', 'MissingPlacement', 'last');
       
    fprintf('Correlation between prev_pr and comp_region\n');
    disp([{'r'} {'r05'} {'r95'} {'pFDR'} {'comp'}; ...
        num2cell(cr(J)) num2cell(cr_lo(J))  num2cell(cr_up(J)) ...
        num2cell(cp(J)) feature_description(J)]);
    

    J = J(cp(J) < 0.05);
    J = flipud(J);

    res = [feature_description(J) num2cell(cr(J)) num2cell(cp(J))]';
    res = sprintf('"%s" (β = %.3f, p = %.3f), ', res{:});
    res = strrep(res, 'p = 0.000', 'p < 0.001');
    res = lower(res);
    disp(res);
    
end
%% Perform comparative analyses between symptom scores
% Compare the observed correlations between the different symptoms (i.e.
% symptom scores). This is done by comparing whether the difference of
% their estimated distribution (obtained using bootstrapping) is
% significantly different from zero.

pres = nan(3,3,nfeatures);

for ifeat = 1:nfeatures

    x = crbtall(ifeat, :, 1) - crbtall(ifeat, :, 2);
    pres(1,2,ifeat) = 2*normcdf(-abs(mean(x) ./ std(x)));
    
    
    x = crbtall(ifeat, :, 3) - crbtall(ifeat, :, 2);
    pres(2,3,ifeat) = 2*normcdf(-abs(mean(x) ./ std(x)));
    
    
    x = crbtall(ifeat, :, 1) - crbtall(ifeat, :, 3);
    pres(1,3,ifeat) = 2*normcdf(-abs(mean(x) ./ std(x)));
    
end

[~, ~, pres_corr] = fdr(pres(:));
pres_corr = reshape(pres_corr, size(pres));

% Report results
fprintf('Corrected across %i tests.\n', nnz(~isnan(pres(:))));
fprintf('Significant differences between features:\n');

for ifeat = 1:nfeatures
    fprintf('%s:\n', feature_description{ifeat});
    for is = 1:nsymptoms
        for iis = 1:nsymptoms
            if pres_corr(is, iis, ifeat) < 0.05

                fprintf(['\tSIGNIFICANT: %s (%.3f) - %s (%.3f) ' ...
                    'p = %.3f (uncorrected p = %.5f)\n'], ...
                    symptoms{is}, mean(crbtall(ifeat, :, is),2), ...
                    symptoms{iis}, mean(crbtall(ifeat, :, iis),2), ...
                    pres_corr(is,iis, ifeat), pres(is,iis, ifeat));
                
            end
        end
    end
end

%% Save all statistics
save(fullfile(results_path, 'results_nke'), 'crall', 'cpall', ...
    'cpall_uncorrected', 'feature_description', 'symptoms', 'pres_corr');

%% Figure presenting results

figure('color', 'white');
colormap(cm);
hold on;

nc = size(crall, 1);
u = repmat([1:nc]', 1, nsymptoms);
v = repmat(1:nsymptoms, nc, 1);

clear h
for i = 1:numel(crall)
    y = u(i) + [-0.5 -0.5 0.5 0.5];
    x = v(i) + [-0.5 0.5 0.5 -0.5];
    h(i) = patch(x, y, crall(i), ...
        'EdgeColor', [1 1 1], ...
        'LineWidth', 0.5);
end

for i = 1:length(h)
    if cpall(i) >= 0.05
        hatch(h(i))
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

axis tight

xticks(1:nsymptoms);
xticklabels(symptoms);
xtickangle(45);
yticks(1:nc);
yticklabels(feature_description);

title('Correlation symptom map with NKE terms')
caxis([-caxis_lim caxis_lim]);

c = colorbar('SouthOutside');
c.Box = false;
c.TickLength = 0;
c.Label.String = 'Correlation (r)';

axis equal

pause(1)
print(fullfile(figures_path, ['nke_grid' resType]), '-dsvg');

diary('off');