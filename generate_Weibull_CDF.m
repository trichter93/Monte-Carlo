function v = generate_Weibull_CDF(k, lambda, I)
    F_v = wblcdf(k, lambda);
    v = generate_X(F_v, I);
end