function X = generate_X_from_inv_CDF(F_Xinv, F_X, I)
    U = rand;
    X = inv_cdf_x_given_I2(U, F_Xinv, F_X, I);
end