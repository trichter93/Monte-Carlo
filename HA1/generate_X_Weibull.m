function X = generate_X_Weibull(lambda, k)
    F_vinv = @(P) wblinv(P, lambda, k);
    X = generate_X(F_vinv);
end