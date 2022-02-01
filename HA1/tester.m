%%

F_X = @(x) normcdf(x);
F_Xinv = @(x) norminv(x);
I=[0.1,2];
Fx
X=zeros(1,1000);
for i = 1 : 1000
    X(i) = generate_X_trunc(F_Xinv, F_X, I);
end

histogram(X);
result = sum(isnan(X(:)));

%% Reality check for Weibull distribution


