function v = generate_X_trunc_Weibull(lambda, k, I)
    F_vinv = @(P) wblinv(P, lambda, k);
    F_v = @(v) wblcdf(v, lambda, k);
    v = generate_X_from_inv_CDF(F_vinv, F_v, I);
end