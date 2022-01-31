%% Beginning of HA1 in Monte Carlo
% For Weibull Distributions, k = Beta = B = shape parameter
% lambda = eta = A = scale parameter
clear;
load powercurve_D236.mat
%% plotting power curve
x = linspace(0,35,1e4);
plot(x, P(x));
grid on;
xlabel('Wind speed v [m/s]')
ylabel('Power output P (v) [kW]');
hold on;

%% constants
rng(0);
conf95 = 1.96;  % 95% confidence interval
d = 236; % height of wind plant
rho = 1.225; % air density
ideal_power_coeff = 16/27; % self explanatory
lambda = [11.7, 10.7, 10.1, 8.8, 8.6, 8.9, 8.6, 8.9, 10.0, 10.9, 11.7, 11.7]; %  Weibull parameters according to HA1 description
k = [2.0, 2.0, 2.0, 1.9, 1.9, 1.9, 1.9, 1.9, 2.0, 1.9, 2.0, 2.0]; %  Weibull parameters according to HA1 description
samples_per_month = 33*[31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]; % Assuming 1 sample per day, can easily be increased
total_N = sum(samples_per_month); % Total number of samples (obviously)

%% Sampling from non-truncated distribution and using CLT to determine Tau_N for each month(amount of power generated) = tau_N, solves first part of a)
V = zeros(12, 31); % Pre-allocating for wind speed samples, row is month, column is day
power_outputs = zeros(12, 31); % Pre-allocating for power output samples, row is month, column is day

% Instead of generate_from_Weibull_inv_cdf use non truncated version:




%% Sampling from truncated Weibull distribution and using CLT to determine Tau_N for each month(amount of power generated) = tau_N, solves second part of a)
V_trunc=zeros(12, 31); % Pre-allocating for wind speed samples, row is month, column is day
power_outputs_trunc = zeros(12, 31); % Pre-allocating for power output samples, row is month, column is day
tau_N_truncated=zeros(1,12); % Pre-allocating for mean power outputs per month
sigma2_truncated=zeros(1,12);
intervals_truncated = zeros(2,12);
I=[3,30]; % the truncated interval
for i = 1 : 12
    for j = 1 : samples_per_month(i)
        V_trunc(i,j) = generate_from_Weibull_inv_CDF(lambda(i), k(i), I);
        power_outputs_trunc(i,j) = P(V_trunc(i,j));
    end
    power_outputs_i = power_outputs_trunc(i,:); % Extracting row of power outputs for each month (Done in order to ignore zeros in next line)
    tau_N_truncated(i) = mean(power_outputs_i(power_outputs_i~=0)); %Taking the mean for each month, ignoring zeros
    sigma2_truncated(i) = var(power_outputs_i(power_outputs_i~=0));
    intervals_truncated(:,i)=[tau_N_truncated(i) - conf95*sqrt(sigma2_truncated(i)/samples_per_month(i)); tau_N_truncated(i) + conf95*sqrt(sigma2_truncated(i)/samples_per_month(i))];
end

tau_N_total = mean(tau_N_truncated);
sigma2_total = var(power_outputs_trunc(:));

interval_total = [tau_N_total - conf95*sqrt(sigma2_total/total_N), tau_N_total + conf95*sqrt(sigma2_total/total_N)];
%%

Ptot = @(v) 1/2 * rho * pi * d.^2 /4 * v.^3; % The analytical function represented by

%% tests

F_X = @(x) normcdf(x);

f_X = @(x) normpdf(x);

I=[0.25,0.75];

Ptest = cdf_x_given_I(0.7, F_X, I);

p= @(x) pdf_x_given_I(x, f_X, I);

x=linspace(0.25,0.75,1e3);

plot(x,p(x));

total = integral(p, 0.27, 0.74,'ArrayValued', true);
%% Clearing all variables except the relevant ones for the rest of the experiment

%clearvars -except [conf95, k, lambda, tau_N];


%% The cdf of X ≤ x given that X in interval I
function F = cdf_x_given_I(x, F_X, I) % 1
    if (x >= I(1) && x<= I(2))
        F = (F_X(x)-F_X(I(1)))/(F_X(I(2))-F_X(I(1)));
    else
        F = 0;
    end
end
%% The pdf of X ≤ x given that X in interval I
function f = pdf_x_given_I(x, f_X, I) % 2
    if(I(1)<=x & x <= I(2))
        norm_constant = integral(f_X, I(1), I(2));
        f=f_X(x)./norm_constant;
    else
        f=0;
    end
end

%% The inverse cdf of X ≤ x given that X in interval I
% If p=0.2 and I=[0.25,0.75], return the x value such that F = cdf_x_given_I(x, F_X, I) = 0.2
function Finv = inv_CDF_x_given_I(U, F_X, I) %Depreciated % 3
    if(U==0)
        Finv = I(1);
    elseif (U==1)
        Finv = I(2);
    else
        points=1e2;
        x = linspace(I(1), I(2), points);
        u = linspace(0, 1, points);
        dp=1e-3;
        for i=1:points
            if(U-cdf_x_given_I(x(i), F_X, I)<=dp)
                %Finv = (u(i) - U)*x(i-1)./(u(i)-u(i-1)) + (U - u(i-1))*x(i)./(u(i)-u(i-1));
                Finv = x(i);
                return;
            end
        end
    end
end

function X = generate_X_from_CDF(F_X,I) %Depreciated 4
    U = rand;
    X = inv_CDF_x_given_I(U, F_X, I);
end

function v = generate_from_Weibull_CDF(k, lambda, I) %Depreciated 5
    F_v = @(v) wblcdf(v, k, lambda);
    v = generate_X_from_inv_CDF(F_v, I);
end

function v = generate_from_Weibull_inv_CDF(lambda, k, I) % 6
    F_vinv = @(P) wblinv(P, lambda, k);
    F_v = @(v) wblcdf(v, lambda, k);
    v = generate_X_from_inv_CDF(F_vinv, F_v, I);
end


function Finv = inv_CDF_x_given_I2(U, F_Xinv, F_X, I) % 7
    if(U==0)
        Finv = I(1);
    elseif (U==1)
        Finv = I(2);
    else
        %Finv = ; %rescaled and so that we exploit curvature of F_Xinv
        Finv = F_Xinv(U*(F_X(I(2)) -  F_X(I(1))) + F_X(I(1)));
    end
end

function X = generate_X_from_inv_CDF(F_Xinv, F_X, I) % 8
    U = rand;
    X = inv_CDF_x_given_I2(U, F_Xinv, F_X, I);
end

