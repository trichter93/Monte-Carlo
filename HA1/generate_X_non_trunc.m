function x = generate_X_non_trunc(FXinv)
    U = rand;
    x = FXinv(U);
end