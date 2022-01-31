function v = generate_from_Weibull_CDF(k, lambda, I)
    F_v = @(v) wblcdf(v, k, lambda);
    v = generate_X_from_CDF(F_v, I);
end