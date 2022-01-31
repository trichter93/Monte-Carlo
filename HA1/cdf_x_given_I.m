function P = CDF_x_given_I(x, F_X, I)
    if (x >= I(1) & x<= I(2))
        P = (F_X(x)-F_X(I(1)))/(F_X(I(2))-F_X(I(1)));
    else
        P = 0;
    end
end

