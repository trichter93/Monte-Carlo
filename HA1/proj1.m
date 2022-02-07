%% Beginning of HA1 in Monte Carlo
% For Weibull Distributions, k = Beta = B = shape parameter
% lambda = eta = A = scale parameter
%clear;
load powercurve_D236.mat
%% plotting power curve
close all;
x = linspace(0,35,1e4);
plot(x, P(x));
grid on;
xlabel('Wind speed v [m/s]')
ylabel('Power output P (v) [W]');
hold on;

%% Power production of a turbine

%constants

rng(0);
conf95 = 1.96;  % 95% confidence interval
d = 236; % height of wind plant
rho = 1.225; % air density
ideal_power_coeff = 16/27; % self explanatory
lambda = [11.7, 10.7, 10.1, 8.8, 8.6, 8.9, 8.6, 8.9, 10.0, 10.9, 11.7, 11.7]; %  Weibull parameters according to HA1 description
k = [2.0, 2.0, 2.0, 1.9, 1.9, 1.9, 1.9, 1.9, 2.0, 1.9, 2.0, 2.0]; %  Weibull parameters according to HA1 description
samples_per_day=33;
samples_per_month = samples_per_day*[31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]; % Assuming some nr of samples per day, can easily be increased
total_N = sum(samples_per_month); % Total number of samples (obviously)

%% Sampling from non-trunc distribution and using CLT to determine Tau_N for each month(amount of power generated) = tau_N, solves first part of a)
V = zeros(12, 31*samples_per_day); % Pre-allocating for wind speed samples, row is month, column is day
power_outputs = zeros(12, 31*samples_per_day); % Pre-allocating for power output samples, row is month, column is day

% Instead of generate_X_trunc_Weibull use non trunc version:

tau_N = zeros(1,12);% Pre-allocating for mean power outputs per month
V_N = zeros(1,12);
sigma2 = zeros(1,12);
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

interval_total = [tau_N_total - conf95*sqrt(sigma2_total/total_N); tau_N_total + conf95*sqrt(sigma2_total/total_N)];



%% Sampling from trunc Weibull distribution and using CLT to determine Tau_N for each month(amount of power generated) = tau_N, solves second part of a)
V_trunc=zeros(12, 31*samples_per_day); % Pre-allocating for wind speed samples, row is month, column is day
power_outputs_trunc = zeros(12, 31*samples_per_day); % Pre-allocating for power output samples, row is month, column is day
tau_N_trunc=zeros(1,12); % Pre-allocating for mean power outputs per month
V_N_trunc = zeros(1,12);
sigma2_trunc=zeros(1,12);
intervals_trunc = zeros(2,12);
I=[3,30]; % the trunc interval
for i = 1 : 12
    for j = 1 : samples_per_month(i)
        V_trunc(i,j) = generate_X_trunc_Weibull(lambda(i), k(i), I);
        power_outputs_trunc(i,j) = P(V_trunc(i,j));
    end % Extracting row of power outputs for each month (Done in order to ignore zeros in next line)
    tau_N_trunc(i) = mean(power_outputs_trunc(i, 1:samples_per_month(i))); %Taking the mean for each month, ignoring zeros
    V_N_trunc(i) = mean(V_trunc(i, 1:samples_per_month(i)));
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
sigma2_CV = zeros(12,1); % fomula : sigma2 +2 * lambda_CV * cov + lambda_CV^2 + var(Y)
tau_CV = zeros(12,1);
lambda_CV=zeros(12,1);
intervals_CV = zeros(2,12);
for i=1:12
    V_mean_exact(i) = gamma(1+1./k(i))*lambda(i);
    covsi = cov(power_outputs_trunc(i, 1:samples_per_month(i)), V_trunc(i,1:samples_per_month(i)));
    covs(i) = covsi(2,1);
    vars(i) = var(V_trunc(i,1:samples_per_month(i)));
    lambda_CV(i) = - covs(i)/vars(i);
    sigma2_CV(i) = sigma2(i) + 2.* lambda_CV(i) * covs(i) + lambda_CV(i).^2 .* var(V_trunc(i,1:samples_per_month(i)));
    tau_CV(i) = tau_N_trunc(i) + lambda_CV(i) * (V_N_trunc(i) - V_mean_exact(i));
    intervals_CV(:,i) = [tau_CV(i) - conf95*sqrt(sigma2_CV(i)/samples_per_month(i)); tau_CV(i) + conf95*sqrt(sigma2_CV(i)/samples_per_month(i))];
end




%% Using importance sampling, ie

V_IS = zeros(12, 31*samples_per_day); % Pre-allocating for wind speed samples, row is month, column is day
power_outputs_IS = zeros(12, 31*samples_per_day); % Pre-allocating for power output samples, row is month, column is day

tau_N_IS=zeros(1,12);% Pre-allocating for mean power outputs per month
V_IS_N = zeros(1,12);
sigma2_IS=zeros(1,12);
intervals_IS = zeros(2,12);
V_var_exact = zeros(1,12);
weights_IS=zeros(12,samples_per_day*31);
vx=linspace(0,35,1000);
for i = 1 : 12
    V_var_exact(i)= wblvar(lambda(i),k(i));
    for j = 1 : samples_per_month(i)
        wbl = @(v) wblpdf(v, lambda(i), k(i));
        [maxy, maxx] = max(P(vx)'.*wbl(vx));
        normtemp = @(v) normcdf(v, vx(maxx), sqrt(V_var_exact(i)));
        norminvtemp = @(v) norminv(v, vx(maxx), sqrt(V_var_exact(i)));
        V_IS(i, j) = generate_X_trunc(norminvtemp, normtemp, I);
        weights_IS(i,j) = wblpdf(V_IS(i,j),lambda(i), k(i))./normpdf(V_IS(i,j), vx(maxx), sqrt(V_var_exact(i)));
        power_outputs_IS(i,j) = P(V_IS(i,j));
    end
    tau_N_IS(i) = mean(power_outputs_IS(i, 1:samples_per_month(i)).*weights_IS(i,1:samples_per_month(i))); %Taking the mean for each month, ignoring zeros
    V_IS_N(i) = mean(V_IS(i, 1:samples_per_month(i)));
    sigma2_IS(i) = var(power_outputs_IS(i, 1:samples_per_month(i)));
    intervals_IS(:,i)=[tau_N_IS(i) - conf95*sqrt(sigma2_IS(i)/samples_per_month(i)); tau_N_IS(i) + conf95*sqrt(sigma2_IS(i)/samples_per_month(i))];
end

tau_N_IS_total = mean(tau_N_IS);
sigma2_IS_total = var(power_outputs_IS(:));

interval_IS_total = [tau_N_IS_total - conf95 * sqrt(sigma2_IS_total/total_N), tau_N_IS_total + conf95*sqrt(sigma2_IS_total/total_N)];


%% Doing Antithetic sampling, ie
%% Doing Antithetic sampling

intervals_AS = zeros(2,12);
var_AS = zeros(1,12);
mu_AS = zeros(1,12);

for i = 1:12
    [var_AS(i), mu_AS(i)] = antithetic(P, @(p) X_trunc(p, @(p1) wblinv(p1, lambda(i), k(i)), @(x1) wblcdf(x1, lambda(i), k(i)), I), samples_per_month(i));
    intervals_AS(:,i) = [mu_AS(i) - conf95*sqrt(var_AS(i)/samples_per_month(i)); mu_AS(i) + conf95*sqrt(var_AS(i)/samples_per_month(i))];
end
%% Calculating probability that turbine delivers power, ie P(power > 0)

PPg0 = zeros(1,12);

for i = 1 : 12
    temp_power_outputs =power_outputs(i, 1:samples_per_month(i));
    PPg0(i) = length(temp_power_outputs(temp_power_outputs ~= 0))/length(temp_power_outputs);
end
%% Estimate average ratio of actual wind turbine output to total wind power, ie E(p)/E(ptot)

Ptot = @(v) 1/2 * rho * pi * d.^2 /4 * v.^3; % The analytical function represented by

% (Use the analytical function Ptot)

average_power_coefficient = zeros(1,12);
intervals_apc = zeros(2,12);
for i =1:12
    average_power_coefficient(i) = tau_N(i)/Ptot(V_N(i));
    intervals_apc(:,i) = [average_power_coefficient(i) - conf95*sqrt(sigma2(i)/samples_per_month(i))/Ptot(V_N(i)), average_power_coefficient(i) + conf95*sqrt(sigma2(i)/samples_per_month(i))/Ptot(V_N(i))]';
end



%% Capacity factor and availability factor

cf = zeros(1,12);

for i = 1 : 12
    cf(i) = sum(power_outputs(i,1:samples_per_month(i)))/(15e6*samples_per_month(i));
end

%% 3 Combined power production of two wind turbines
%% tests

% F_X = @(x) normcdf(x);
% 
% f_X = @(x) normpdf(x);
% 
% I=[0.25,0.75];
% 
% Ptest = F_X_trunc(0.7, F_X, I);
% 
% p= @(x) fX_trunc(x, f_X, I);
% 
% x=linspace(0.25,0.75,1e3);
% 
% plot(x,p(x));
% 
% total = integral(p, 0.27, 0.74,'ArrayValued', true);
%% Clearing all variables except the relevant ones for the rest of the experiment

%clearvars -except [conf95, k, lambda, tau_N];

%% Exercise 3, Multivariate Weibull distributions
%constants
alpha3 = 0.638;
p3 = 3;
q3 = 1.5;

lambda3 = 10.05;
k3 = 1.95;

%% Defining the pdf, cdf, bivariate pdf and bivariate cdf

Fv = @(v) wblcdf(v, lambda3, k3);
fv = @(v) wblpdf(v, lambda3, k3);


F12 = @(v1,v2) Fv(v1) .* Fv(v2) .* (1 + alpha3 .* (1-Fv(v1).^p3).^q3 .* (1-Fv(v2).^p3).^q3);
f12 = @(v1,v2)  fv(v1) .* fv(v2) .* (1 + alpha3.*(1-Fv(v1).^p3).^(q3-1) .* (1-Fv(v2).^p3).^(q3-1) .* (Fv(v1).^p3 .* (1+p3.*q3) - 1) .* (Fv(v2).^p3 .* (1+p3.*q3) - 1));

%integral2(f12, 0, 10, 0, 10) % sanity check
%F12(10,10) % should and is equal to one line above
%% a, ie finding E(P(V1)+P(V2))
% Each marginal distribution is Weibull, and E(P(V1) = E(P(V2)). Hence only
% 1 needs to be calculated. => Use importance sampling again

V3 = zeros(1, 1000); % Pre-allocating for wind speed samples, row is month, column is day
power_outputs3 = zeros(1, 1000); % Pre-allocating for power output samples, row is month, column is day

% Instead of generate_X_trunc_Weibull use non trunc version:
% Pre-allocating for mean power outputs per month

weights3=zeros(1,1000);
vx=linspace(0,35,1000);
V_var_exact3 = wblvar(lambda3,k3);
wbl3 = @(v) wblpdf(v, lambda3, k3);
[maxy, maxx] = max(P(vx)'.*wbl3(vx));
normtemp = @(v) normcdf(v, vx(maxx), sqrt(V_var_exact3));
norminvtemp = @(v) norminv(v, vx(maxx), sqrt(V_var_exact3));
for j = 1 : 1000
    V3(j) = generate_X_trunc(norminvtemp, normtemp, I);
    weights3(j) = wblpdf(V3(j), lambda3, k3)./normpdf(V3(j), vx(maxx), sqrt(V_var_exact3));
    power_outputs3(j) = P(V3(j));
end
tau_N3 = mean(power_outputs3(:).*weights3(:));
V_N3 = mean(V3(:));
sigma23 = var(power_outputs3(:));
interval3=[tau_N3 - conf95*sqrt(sigma23/1000); tau_N3 + conf95*sqrt(sigma23/1000)];
%% b, ie finding C(P(v1), P(v2)). Use classical expression C(X,Y) = E(XY) - E(X)E(Y)

% yields 3 integrals:

% E(X) = integral( x * f(x,y) dx dy), E(Y)=integral( y * f(x,y) dx dy) KNOWN
% E(XY) = integral( xy * f(x,y) dx dy)

% Should use importance sampling for both directions

% XY = P(V1)*P(V2) where V1 and V2 are generated with the help of mvnrnd 

% In order to generate a 2D space that can be used to find the values
% P(v1,v2)*wblpdf(v1,v2) Use meshgrid on two of the same size linspaces
V32 = zeros(2, 1000); % Pre-allocating for wind speed samples
power_outputs32 = zeros(2, 1000); % Pre-allocating for power output samples


% where mu is argmax( P

v1 = linspace(0,35,1000);
v2 = linspace(0,35,1000);

[V1,V2] = meshgrid(v1,v2);

weights32 = zeros(1,1000);
f12dist = f12(V1, V2);
PV1 = reshape(P(V1),[1000,1000]);
PV2 = reshape(P(V2),[1000,1000]);
prod_mat = ((PV1'.*PV2)'.*f12dist)'; % the values per row and the columns they occur at
[~, row] = max(max(prod_mat));
[~, col] = max(prod_mat(row,:));
covar_mat = 20*eye(2);
for i = 1:1000
    V32(:,i)=mvnrnd(mu, covar_mat);
    weights32(i) = f12(V32(1,i),V32(2,i))/mvnpdf([V32(1,i),V32(2,i)], mu, covar_mat);
    power_outputs32(:,i)=P(V32(:,i));
end

E_PV1PV2 = mean(power_outputs32(1,:) .* power_outputs32(2,:).*weights32);

C_PV1PV2 = E_PV1PV2 - tau_N3^2;
%% The variance and the std of the sum P(V1)+P(V2), perhaps too simple... But no idea why it wouldnt be like this

Var_PV1pPV2 = 2 * C_PV1PV2 + 2 * sigma23;

std_PV1pPV2 = sqrt(Var_PV1pPV2);

%% Finding the 95% conf interval for the probabilities P(P(V1)+P(V2) > 15 MW) and P(P(V1)+P(V2) < 15 MW)

% If I can sample from this distribution then that would surely help

% I have F(v1,v2) and f(v1,v2), ie P(V1 < v1 and V2 < v2) and the corresponding pdf
N= 10000;
V3d2 = zeros(2, N); % Pre-allocating for wind speed samples
power_outputs3d2 = zeros(2, N); % Pre-allocating for power output samples
weights3d2 = zeros(1,N);
successes = zeros(1,N);
fails = ones(1,N);
for i = 1:N
    V3d2(:,i)=mvnrnd(mu, covar_mat);
    weights3d2(i) = f12(V3d2(1,i),V3d2(2,i))/mvnpdf([V3d2(1,i),V3d2(2,i)], mu, covar_mat);
    power_outputs3d2(:,i)=P(V3d2(:,i));
    if sum(power_outputs3d2(:,i)) > 15*1e6
        successes(i) = 1;
        fails(i) = 0;
    end
end
var_3d2s = var(successes.*weights3d2);
P_3d2s = mean(successes.*weights3d2);

var_3d2f = var(fails.*weights3d2);
P_3d2f = mean(fails.*weights3d2);


interval_3d2s = [P_3d2s - conf95 * sqrt(var_3d2s/N) ; P_3d2s + conf95 * sqrt(var_3d2s/N)];

interval_3d2f = [P_3d2f - conf95 * sqrt(var_3d2f/N) ; P_3d2f + conf95 * sqrt(var_3d2f/N)];

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

function m = wblvar(lambda, k)
    m = lambda.^2 * (gamma(1+2/k)-gamma(1+1/k).^2);
end
function [vAS, muAS] = antithetic(phi, Finv, N)
    U = rand(N,1); %uniform distribustion

    V = phi(Finv(U));
    Vt = phi(Finv(1-U));

    mu1 = (1/N)*sum(V); 
    mu2 = (1/N)*sum(Vt);

    muAS = (mu1+mu2)/2;

    covM = cov(V,Vt);
    vAS = sum(covM(:))/4;
end