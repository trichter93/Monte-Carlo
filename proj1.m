%% Beginning of HA1 in Monte Carlo

load powercurve_D236.mat
%%
x = linspace(0,35,1e4);

plot(x, P(x));
grid on;
xlabel('Wind speed v [m/s]')
ylabel('Power output P (v) [kW]');
hold on;
%% 
d = 236; 
rho = 1.225;
ideal_power_coeff = 16/27;



%%

Ptot = @(v) 1/2 * rho * pi * d.^2 /4 * v.^3;

%% tests

F_X = @(x) normcdf(x);

f_x = @(x) normpdf(x);

I=[0.25,0.75];

Ptest = cdf_x_given_I(0.7, F_x, I);

p= @(x) pdf_x_given_I(x, f_x, I);

x=linspace(0.25,0.75,1e3);

plot(x,p(x));

total = integral(p, 0.27, 0.74,'ArrayValued', true);

inv_cdf_x_given_I(0.3, F_X, I)



%% The cdf of X ≤ x given that X in interval I
function F = cdf_x_given_I(x, F_X, I)
    if (x >= I(1) && x<= I(2))
        F = (F_X(x)-F_X(I(1)))/(F_X(I(2))-F_X(I(1)));
    else
        F = 0;
    end
end
%% The pdf of X ≤ x given that X in interval I
function f = pdf_x_given_I(x, f_X, I)
    if(I(1)<=x & x <= I(2))
        norm_constant = integral(f_X, I(1), I(2));
        f=f_X(x)./norm_constant;
    else
        f=0;
    end
end

%% The inverse cdf of X ≤ x given that X in interval I
% If p=0.2 and I=[0.25,0.75], return the x value such that F = cdf_x_given_I(x, F_X, I) = 0.2
function Finv = inv_cdf_x_given_I(U, F_X, I)
    if(U==0)
        Finv = I(1);
    elseif (U==1)
        Finv = I(2);
    else
        points=1e4;
        x = linspace(I(1), I(2), points);
        u = linspace(0, 1, points);
        dp=1e-5;
        for i=1:points
            if(abs(U-cdf_x_given_I(x(i), F_X, I))<=dp)
                Finv = (u(i+1) - U)*x(i)/(u(i+1)-u(i)) + (U - u(i))*x(i+1)/(u(i+1)-u(i));
                return;
            end
        end
    end
end

function X = generate_X(F_X,I)
    U = rand;
    X = inv_cdf_x_given_I(U, F_X, I);
end

