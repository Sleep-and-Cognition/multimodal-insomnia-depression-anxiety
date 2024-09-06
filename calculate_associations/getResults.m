function [res, res_perm, res_bt] = getResults(this_data, YY, covariates, nperm, nboot, nsample)
%getResults Calculate the symptom maps and null models for the given data.
% INPUTS:
% thisData      - Region/connectivity data (regions x subjects)
% YY            - Symptom scores (subjects x symptoms)
% covariatesAll - Covariate data (subjects x covariates)
% nperm         - Number of permutations for the null model
% nboot         - Number of bootstraps for estimating confidence intervals
% nsample       - Number of subjects to sample for bootstrapping (0 for all)
%
% OUTPUTS:
% res      - Calculated results for each region/connection (regions x symptoms)
% res_perm - Calculated permutation-sample results for each
%            region/connection (regions x symptoms x permutations)
% res_bt   - Calculated bootstrap-sample results for each 
%            region/connection (regions x symptoms x bootstraps)

% NOTE:
% The implementation of this function is validated in test_getResults.m

nsymptoms = size(YY,2);

% Initialize result matrices
res = nan(size(this_data, 1), nsymptoms);
res_perm = nan(size(this_data, 1), nsymptoms, nperm);
res_bt = nan(size(this_data, 1), nsymptoms, nboot);

% Convert data to single precision
XX = single(this_data'); % Predictors
CC = single(covariates); % Covariates
YY = single(YY); % Outcome variables

% Exclude subjects with missing symptom scores (YY), covariates or no
% predictor data (XX) (ie. subjects excluded in QC step).
indx_others = ~any(isnan(YY(:, 1:nsymptoms)),2) & ~any(isnan(covariates),2) ...
    & any(~isnan(XX),2);
XX = XX(indx_others,:);
CC = CC(indx_others,:);
YY = YY(indx_others,:);
fprintf('Max number of subjects included in model: %i\n', nnz(indx_others));

nSubjects = nnz(indx_others);

% limit bootstrap sample for bootstrapping during sensitivty analysis.
if nsample <= 0
    nsample = nSubjects;
end

fprintf('Initialize permutations\n');
tic
% Define permutations of the null model. (Define them in advance such that
% same permutations are used for each disorder and brain measure to ensure
% that their covariances are preserved).
indx_perms = nan(nSubjects, nperm);
for i = 1:nperm
    indx_perms(:,i) = randperm(nSubjects,nSubjects);
end

% Define permutations for bootstrapping of estimate.

indx_boots = randi(nSubjects, nsample, nboot);
toc

fprintf('Run main and permutations\n');
tic

% Loop through each symptom
for isymptom = 1:nsymptoms
        
    Y = YY(:, isymptom);
    
    % Loop through each region/connection
    for ir = 1:size(this_data, 1)
        
        X = XX(:, ir);
        
        % indx selects all subjects in which the brain measure is available
        % (e.g. in DWI a connection might be missing in a percentage of the
        % people).
        indx = ~isnan(X);
        nindx = nnz(indx);

        % Perform regression model.
        y = Y(indx);
        y = y ./ std(y);
        data = [X(indx) CC(indx,:)];
        data = data ./ std(data, [], 1);
        [Q,R] = qr([ones(size(data,1),1) data],0);
        beta = R\(Q'*y);
        res(ir, isymptom) = beta(2);
        
        % Perform permutations for the null model
        % The outcome variable is permuted, therefore the variables defined
        % by the predictors and covariates (indx, Q and R) can be reused,
        % which speeds up the analyses (by a lot).
        for ii = 1:nperm           
            
            this_perm = indx_perms(:,ii);
            
            y = Y(this_perm);
            y = y(1:nindx);
            y = y ./ std(y);
            beta = R\(Q'*y);
            res_perm(ir, isymptom, ii)  = beta(2);    
            
        end
    end
end
toc

%% Perform bootstrapping to estimate confidence intervals.
% For debugging:
% regstats(y(indx_boots(:, ii)), data(indx_boots(:, ii), :), 'linear', {'beta'});

fprintf('Run bootstrap\n');
tic
% Loop through each bootstrap sample
for ii = 1:nboot

    % Loop through each region/connection
    for ir = 1:size(this_data, 1)

        X = XX(:, ir);

        this_boot = indx_boots(:,ii);

        indx = ~isnan(X(this_boot));
        indx = this_boot(indx);
        nindx = nnz(indx);

        data = [X(indx) CC(indx,:)];
        data = data ./ std(data, [], 1);
        [Q,R] = qr([ones(size(data,1),1) data],0);

        % Loop through each symptom
        for isymptom = 1:nsymptoms

            Y = YY(:, isymptom);

            y = Y(indx);
            y = y ./ std(y);

            beta = R\(Q'*y);
            res_bt(ir, isymptom, ii)  = beta(2);

        end

    end
end
toc