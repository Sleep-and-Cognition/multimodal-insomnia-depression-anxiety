%% Replace subject-label permutation with spin-based permutation

%% Initialize
data = load(fullfile(absolute_path, ...
    ['nke_data/data/hansen_receptors/spins-' atlas '.mat']));

[~, I, J] = intersect(regionDescriptions, data.regions, 'stable');
assert(isequal(regionDescriptions(I)', data.regions(J)));

spins = nan(length(regionDescriptions), size(data.spins,2));
spins(I,:) = data.spins(J,:);

%% Compute permutations using spin-based permutation
% Subcortical regions do not have spin-based permutations, instead they are
% being randomly suffled.

% Structural and functional connectivity
if contains(data_modality, {'connvF', 'connvS'})

    spins = spins + 15; % Convert from python indexing

    % Randomly shuffle subcortical regions
    indx_perms = nan(14, size(data.spins,2));
    for i = 1:size(data.spins,2)
        indx_perms(:,i) = randperm(14);
    end

    spins(1:14, :) = indx_perms;

    % Apply shuffle and spin-based permutations
    difference_t_perm = nan(size(difference_t, 1), ...
        nsymptoms, size(spins,2));

    for i = 1:3
        A = squareform(difference_t(:, i));
        for ii = 1:size(spins, 2)
            difference_t_perm(:, i, ii) = squareform(A(spins(:, ii), ...
                spins(:, ii)));
        end
    end

% Subcortical volumes
elseif strcmp(data_modality, 'volumes')

    % Randomly shuffle subcortical regions
    indx_perms = nan(14, size(data.spins,2));
    for i = 1:size(data.spins,2)
        indx_perms(:,i) = randperm(14);
    end

    spins(1:14, :) = indx_perms;

    % Apply shuffle-based permutations
    difference_t_perm = nan(length(regionDescriptions), ...
        nsymptoms, size(spins,2));
    for i = 1:3
        for ii = 1:size(spins,2)
            difference_t_perm(:, i, ii) = difference_t(spins(:,ii), i);
        end
    end

% Cortical surface area and cortical thickness
else

    % Convert spins from python indexing
    spins = spins + 1;

    % Apply spin-based permutations
    difference_t_perm = nan(length(regionDescriptions), ...
        nsymptoms, size(spins,2));

    for i = 1:3
        for ii = 1:size(spins,2)
            difference_t_perm(:, i, ii) = difference_t(spins(:,ii), i);
        end
    end
end