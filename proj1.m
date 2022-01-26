%% Beginning of HA1 in Monte Carlo

load powercurve_D236.mat  

%% 


%% tests

F_X = @(x) unifcdf(x);

f_x = @(x) unifpdf(x);

P = cdf_x_given_I(0.7, F_x, [0.25,0.75]);

p= @(x) pdf_x_given_I(x, f_x, [0.25,0.75]);

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
function Finv = inv_cdf_x_given_I(p, F_X, I)
    if(p==0)
        Finv = I(1);
    elseif (p==1)
        Finv = I(2);
    else
        points=1e4;
        x = linspace(I(1), I(2), points);
        dp=1e-5;
        for i=1:points
            if(p-cdf_x_given_I(x(i), F_X, I)<=dp)
                Finv = x(i);
                return;
            end
        end
    end
end
