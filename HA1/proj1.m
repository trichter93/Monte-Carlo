%% Beginning of HA1 in Monte Carlo
% For Weibull Distributions, k = Beta = B = shape parameter
% lambda = eta = A = scale parameter
%clear;
load powercurve_D236.mat
%% plotting power curve
x = linspace(0,35,1e4);
plot(x, P(x));
grid on;
xlabel('Wind speed v [m/s]')
ylabel('Power output P (v) [kW]');
hold on;

%% Power production of a turbine

%constants

rng(0);
conf95 = 1.96;  % 95% confidence interval
d = 236; % height of wind plant
rho = 1.225; % air density
ideal_power_coeff = 16/27; % self explanatory
samples_per_day=33;
lambda = [11.7, 10.7, 10.1, 8.8, 8.6, 8.9, 8.6, 8.9, 10.0, 10.9, 11.7, 11.7]; %  Weibull parameters according to HA1 description
k = [2.0, 2.0, 2.0, 1.9, 1.9, 1.9, 1.9, 1.9, 2.0, 1.9, 2.0, 2.0]; %  Weibull parameters according to HA1 description
samples_per_month = samples_per_day*[31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]; % Assuming 1 sample per day, can easily be increased
total_N = sum(samples_per_month); % Total number of samples (obviously)

%% Sampling from non-trunc distribution and using CLT to determine Tau_N for each month(amount of power generated) = tau_N, solves first part of a)
V = zeros(12, 31*samples_per_day); % Pre-allocating for wind speed samples, row is month, column is day
power_outputs = zeros(12, 31*samples_per_day); % Pre-allocating for power output samples, row is month, column is day

% Instead of generate_X_trunc_Weibull use non trunc version:

tau_N=zeros(1,12);% Pre-allocating for mean power outputs per month
V_N = zeros(1,12);
sigma2=zeros(1,12);
intervals = zeros(2,12);
for i = 1 : 12
    for j = 1 : samples_per_month(i)
        V(i,j) = wblrnd(lambda(i), k(i));
        power_outputs(i,j) = P(V(i,j));
    end
    tau_N(i) = mean(power_outputs(i, 1:samples_per_month(i))); %Taking the mean for each month, ignoring zeros
    V_N(i) = mean(V(i, 1:samples_per_month(i)));
    sigma2(i) = var(power_outputs(i, 1:samples_per_month(i)));
    intervals(:,i)=[tau_N(i) - conf95*sqrt(sigma2(i)/samples_per_month(i)); tau_N(i) + conf95*sqrt(sigma2(i)/samples_per_month(i))];
end

tau_N_total = mean(tau_N);
sigma2_total = var(power_outputs(:));

interval_total = [tau_N_total - conf95*sqrt(sigma2_total/total_N), tau_N_total + conf95*sqrt(sigma2_total/total_N)];



%% Sampling from trunc Weibull distribution and using CLT to determine Tau_N for each month(amount of power generated) = tau_N, solves second part of a)
V_trunc=zeros(12, 31*samples_per_day); % Pre-allocating for wind speed samples, row is month, column is day
power_outputs_trunc = zeros(12, 31*samples_per_day); % Pre-allocating for power output samples, row is month, column is day
tau_N_trunc=zeros(1,12); % Pre-allocating for mean power outputs per month
sigma2_trunc=zeros(1,12);
intervals_trunc = zeros(2,12);
I=[3,30]; % the trunc interval
for i = 1 : 12
    for j = 1 : samples_per_month(i)
        V_trunc(i,j) = generate_X_trunc_Weibull(lambda(i), k(i), I);
        power_outputs_trunc(i,j) = P(V_trunc(i,j));
    end % Extracting row of power outputs for each month (Done in order to ignore zeros in next line)
    tau_N_trunc(i) = mean(power_outputs_trunc(i, 1:samples_per_month(i))); %Taking the mean for each month, ignoring zeros
    sigma2_trunc(i) = var(power_outputs_trunc(i, 1:samples_per_month(i)));
    intervals_trunc(:,i)=[tau_N_trunc(i) - conf95*sqrt(sigma2_trunc(i)/samples_per_month(i)); tau_N_trunc(i) + conf95*sqrt(sigma2_trunc(i)/samples_per_month(i))];
end

tau_N_total_trunc = mean(tau_N_trunc);
sigma2_total_trunc = var(power_outputs_trunc(:));

interval_total_trunc = [tau_N_total_trunc - conf95*sqrt(sigma2_total_trunc/total_N), tau_N_total_trunc + conf95*sqrt(sigma2_total_trunc/total_N)];
%% Using V as control variate, ie mu^hat_CV = mu^hat_MC + lambda*(theta^hat_MC - theta)
% Everything should still be done for each month
% V_mean (theta) is analytically average wind speed, which is given by Gamma(1+1/k)*lambda

% power_outputs (X_i) corresponds to tau_N (mu^hat_MC) and are the power outputs

% V (Y_i) corresponds to V_N (theta^hat_MC) and are the wind speeds

% find var(tau_N) (= sigma2/n ?)
% find covar(tau_N, V_N) = 1/n^2 * sum(diag(cov(X , Y))) ?
V_mean_exact = zeros(12,1); % The exact mean wind speeds (theta)
vars = zeros(12, 1); % the variances for each month
covs = zeros(12, 1); % The covariances we are interested in for determining the optimal lambda_CV
tau_CV = zeros(12,1);
lambda_CV=zeros(12,1);
for i=1:12
    V_mean_exact(i) = gamma(1+1./k(i))*lambda(i);
    covs(i) = my_cov(power_outputs(i, 1:samples_per_month(i)), V(i,1:samples_per_month(i)));
    vars(i) = my_var(V(i,1:samples_per_month(i)));
    lambda_CV(i) = - covs(i)/vars(i);
    tau_CV(i) = tau_N(i) + lambda_CV(i) * (V_N(i) - V_mean_exact(i));
end



%% Using importance sampling, ie

%% Doing Antithetic sampling, ie

%% Calculating probability that turbine delivers power, ie P(power > 0)

%% Estimate average ratio of actual wind turbine output to total wind power, ie E(p)/E(ptot)

Ptot = @(v) 1/2 * rho * pi * d.^2 /4 * v.^3; % The analytical function represented by

% (Use the analytical function Ptot

%% Capacity factor and availability factor


%% 3 Combined power production of two wind turbines
%% tests

F_X = @(x) normcdf(x);

f_X = @(x) normpdf(x);

I=[0.25,0.75];

Ptest = F_X_trunc(0.7, F_X, I);

p= @(x) fX_trunc(x, f_X, I);

x=linspace(0.25,0.75,1e3);

plot(x,p(x));

total = integral(p, 0.27, 0.74,'ArrayValued', true);
%% Clearing all variables except the relevant ones for the rest of the experiment

%clearvars -except [conf95, k, lambda, tau_N];

%%

%% The cdf of X ≤ x given that X in interval I
function F = F_X_trunc(x, F_X, I) % 1
    if (x >= I(1) && x<= I(2))
        F = (F_X(x)-F_X(I(1)))/(F_X(I(2))-F_X(I(1)));
    else
        F = 0;
    end
end
%% The pdf of X given that X in interval I
function f = fX_trunc(x, f_X, I) % 2
    if(I(1)<=x & x <= I(2))
        norm_constant = integral(f_X, I(1), I(2));
        f=f_X(x)./norm_constant;
    else
        f=0;
    end
end

function v = generate_X_trunc_Weibull(lambda, k, I) % 3
    F_vinv = @(P) wblinv(P, lambda, k);
    F_v = @(v) wblcdf(v, lambda, k);
    v = generate_X_trunc(F_vinv, F_v, I);
end


function Finv = X_trunc(U, F_Xinv, F_X, I) % 4
    if(U==0)
        Finv = I(1);
    elseif (U==1)
        Finv = I(2);
    else
        %Finv = ; %rescaled and so that we exploit curvature of F_Xinv
        Finv = F_Xinv(U*(F_X(I(2)) -  F_X(I(1))) + F_X(I(1)));
    end
end

function X = generate_X_trunc(F_Xinv, F_X, I) % 5
    U = rand;
    X = X_trunc(U, F_Xinv, F_X, I);
end

function var = my_var(X)
    N= length(X);
    mX = mean(X);
    var = 0;
    for i = 1:N
        var = var + (X(i) - mX).^2 / N / (N-1);
    end
end

function cov = my_cov(X,Y)
    N= length(X);
    mX = mean(X);
    mY = mean(Y);
    cov=0;
    for i =1:N
        cov = cov + (X(i) - mX) * (Y(i) - mY)/N/(N-1);
    end
end
