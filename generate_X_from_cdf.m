function X = generate_X_from_CDF(F_X, I)
    U = rand;
    X = inv_cdf_x_given_I(U, F_X, I);
end