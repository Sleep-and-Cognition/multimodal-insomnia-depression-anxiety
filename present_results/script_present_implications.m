%% Estimate the clinical impact of the observed association strengths
% In this script we will estimate the clinical impact of the observed
% association between anxiety symptom severity and cortical surface area.

%% Initialize
% Run the initialization script and any necessary initialization steps. Run
% script_calculate_regional_maps_IDA for the anxiety and cortical surface
% modalities. Also, make sure to run getResults until the line where [Q,R]
% is calculated.

% Now, one can replicate the results from the paper, by running the
% following:

% K>> fitlm(data, y) 
% ans = 
% 
% 
% Linear regression model:
%     y ~ 1 + x1 + x2 + x3 + x4 + x5 + x6 + x7
% 
% Estimated Coefficients:
%                    Estimate        SE         tStat       pValue  
%                    _________    _________    _______    __________
% 
%     (Intercept)       1.1856     0.089653     13.225      8.54e-40
%     x1             -0.049933    0.0076563    -6.5218    7.0744e-11
%     x2              -0.18918      0.09799    -1.9307      0.053537
%     x3               0.24187     0.053878     4.4892    7.1804e-06
%     x4               0.13162     0.098354     1.3383       0.18082
%     x5              -0.23413      0.05436    -4.3071    1.6602e-05
%     x6             -0.022392    0.0092335    -2.4251      0.015312
%     x7             -0.016604    0.0062063    -2.6754     0.0074695
% 
% 
% Number of observations: 25763, Error degrees of freedom: 25755
% Root Mean Squared Error: 0.995
% R-squared: 0.00998,  Adjusted R-Squared: 0.00971
% F-statistic vs. constant model: 37.1, p-value = 0

%% Estimate effect size associated with having any symptom
% In this section, we will estimate the effect size associated with having
% any symptom by using a generalized model that fits a binomial logistic
% regression:

fitglm(data, y>0, 'distribution' ,'binomial')
b = glmfit(data, y>0, 'binomial', 'link', 'logit');

fprintf(['People with %.2f%% (one standard deviation) smaller ' ...
    'CSA are %.2f times more likely to have any anxiety symptom.\n'], ...
    100 / mean(data(:,1)), exp(-b(2)));

% The output should be:

% ans = 
% 
% 
% Generalized linear regression model:
%     logit(y) ~ 1 + x1 + x2 + x3 + x4 + x5 + x6 + x7
%     Distribution = Binomial
% 
% Estimated Coefficients:
%                    Estimate        SE        tStat        pValue  
%                    _________    ________    ________    __________
% 
%     (Intercept)      0.60285     0.18718      3.2208     0.0012785
%     x1             -0.096282    0.016007      -6.015    1.7984e-09
%     x2              -0.19808     0.20566    -0.96317       0.33546
%     x3               0.60458     0.11213      5.3917    6.9791e-08
%     x4              0.053282     0.20685     0.25759       0.79672
%     x5              -0.57811     0.11374     -5.0827    3.7217e-07
%     x6             -0.019417    0.019251     -1.0087       0.31314
%     x7             -0.031617    0.012838     -2.4628      0.013785
% 
% 
% 25763 observations, 25755 error degrees of freedom
% Dispersion: 1
% Chi^2-statistic vs. constant model: 299, p-value = 0

% Comparing the top and bottom 10%:
P = prctile(data(:,1), [10 90]);

r_top = mean(YY(data(:,1) > P(2)) == 1);
r_bottom = mean(YY(data(:,1) < P(1)) == 1);

r_diff = 100*((r_top - r_bottom) ./ r_bottom);

fprintf('percentage with any anxiety symptoms among bottom = %.3f\n', r_bottom);
fprintf('percentage with any anxiety symptoms among top = %.3f\n', r_top);
fprintf('percentage increase = %.2f%% \n', r_diff);
