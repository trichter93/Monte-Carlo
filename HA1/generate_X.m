function x = generate_X(FXinv)
    U = rand;
    x = FXinv(U);
end