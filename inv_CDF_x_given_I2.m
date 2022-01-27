function Finv = inv_CDF_x_given_I2(U, F_Xinv, F_X, I)
    if(U==0)
        Finv = I(1);
    elseif (U==1)
        Finv = I(2);
    else
        %Finv = ; %rescaled and so that we exploit curvature of F_Xinv
        Finv = F_Xinv(U*(F_X(I(2)) -  F_X(I(1)))+ F_X(I(1)));
    end
end