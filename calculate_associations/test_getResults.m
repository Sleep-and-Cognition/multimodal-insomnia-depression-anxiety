% TEST_GETRESULTS Test script for getResults
%
% This script asserts that the function getResults provides results similar
% to regstats and fitlm by comparing their results on synthetic data.
%
% Test 1: Testing beta coefficient
% This test creates synthetic data and calculates the beta coefficient
% using both regstats and getResults functions. It compares the results and
% asserts that the difference is within a tolerance.
%
% Test 2: Testing permutations
% This test creates synthetic data and performs permutations using both
% fitlm and getResults functions. It compares the p-values obtained from
% the two methods and calculates the average p-value and power.
%
% Test 3: Testing bootstrapping for confidence intervals
% This test creates synthetic data and calculates the confidence intervals
% using both fitlm and getResults functions. It compares the confidence
% intervals obtained from the two methods and calculates the average
% difference between the observed and estimated values of r.
%
% All tests are expected to pass without any assertion errors.

%% Test beta coefficient wrt regstats
fprintf('Testing beta coefficient...\n');
tic
nperm = 5000;
nboot = 5000;
nsub = 20000;
nsample = 0; % do not limit sample size in bootstrapping procedure

% Create synthetic data
sigma = [1.0000, -0.029;
    -0.029,  1.0000];
R = mvnrnd([10 10], sigma, nsub);
Cr = randn([nsub 1]);
Xr = R(:,1)';
Yr = R(:,2);

% Calculate reference results
X = zscore([Xr' Cr]);
Y = zscore(Yr);
res_ref = regstats(Y, X, 'linear', {'beta'});
res_ref = res_ref.beta(2);

% Run implementation and compare results
[res, ~, ~] = getResults(Xr, Yr, Cr, nperm, nboot, nsample);
assert(abs(res - res_ref) < 0.0001);

fprintf(' finished successfully.\n');
toc
%% Test permutations with respect to fitlm
fprintf('Testing permutations...\n');

tic

nperm = 1000;
nboot = 2;
nsub = 20000;
nsample_test = 20;
sigma = [1.0000, -0.029;
    -0.029,  1.0000];
nsample = 0; % do not limit sample size in bootstrapping procedure


fprintf('\tOverview results per iteration:\n')
fprintf('\tRef\tObs\n');

results = nan(nsample_test,2);
for i = 1:nsample_test

    % Create synthetic data
    R = mvnrnd([10 10], sigma, nsub);
    Cr = randn([nsub 1]);

    Xr = R(:,1)';
    Yr = R(:,2);

    % Calculate reference results
    X = zscore([Xr' Cr]); % acts on first dimension.
    Y = zscore(Yr);
    
    l = fitlm(X, Y);
    p_ref = l.Coefficients.pValue(2);


    % Run implementation and compare results
    [res, res_perm, ~] = getResults(Xr, Yr, Cr, nperm, nboot, nsample);
    [~, p] = ttest2((res), squeeze(res_perm));
    
    results(i,:) = [p_ref, p];
    fprintf('\t%.4f\t%.4f\n', p_ref, p);

end

fprintf('Average p-val:\n');
fprintf('\t%.4f\t%.4f\n', mean(results));

% Compare power of both methods
results_power = mean(results < (0.05 / 80), 1);

fprintf('Average power:\n');
fprintf('\t%.4f\t%.4f\n', results_power);

fprintf('...finished successfully.\n');
toc

%% Test bootstrapping for confidence intervals
fprintf('Testing bootstrapping for confidence intervals...\n');

tic

nperm = 2;
nboot = 1000;
nsub = 2000;
nsample_test = 20;
sigma = [1.0000, -0.029;
    -0.029,  1.0000];
nsample = 0; % do not limit sample size in bootstrapping procedure

fprintf('\tOverview results per iteration:\n')
fprintf('\tRef\t\t\tObs\n');

results = nan(nsample_test,4);
results2 = nan(nsample_test,2);
for i = 1:nsample_test

    % Create synthetic data
    R = mvnrnd([10 10], sigma, nsub);
    Cr = randn([nsub 1]);

    Xr = R(:,1)';
    Yr = R(:,2);

    % Calculate reference results
    X = zscore([Xr' Cr]);
    Y = zscore(Yr);

    l = fitlm(X,Y);
    tmp = coefCI(l);
    ci_ref = tmp(2,:);

    % Run implementation and compare results
    [~, ~, res_bt] = getResults(Xr, Yr, Cr, nperm, nboot, nsample);
    ci = prctile(res_bt,[2.5 97.5]);


    results(i,:) = [ci_ref, ci];
    results2(i,:) = [mean(res_bt) l.Coefficients.Estimate(2)];
    fprintf('\t[%.4f, %.4f]\t[%.4f, %.4f]\n', results(i,:));
end


results2_diff = mean(abs(results2(:,2) - results2(:,1)));
fprintf('average differene between r observed and r estimated = %.4f (%.1f%%)\n', ...
    results2_diff, 100 * results2_diff ./ mean(abs(results2(:,1))));

results = mean([results(:,2) - results(:,1) results(:,4) - results(:,3)],1);
assert(abs(results(1) - results(2)) < 0.01);

fprintf('...finished successfully.\n');

fprintf('\nAll tests passed.\n');

toc

