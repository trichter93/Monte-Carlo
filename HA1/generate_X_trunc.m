function X = generate_X_trunc(F_Xinv, F_X, I)
    U = rand;
    X = X_trunc(U, F_Xinv, F_X, I);
end