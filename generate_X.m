function X = generate_X(F_X,I)
    U = rand;
    X = inv_cdf_x_given_I(U, F_X,I);
end