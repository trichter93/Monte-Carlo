function X = generate_X_from_inv_CDF(F_Xinv, I)
    U = rand;
    X = inv_cdf_x_given_I2(U, F_Xinv, I);
end