function x = inv_CDF_x_given_I2(P, F_Xinv, F_X, I)
    if(P==0)
        x = I(1);
    elseif (P==1)
        x = I(2);
    else
        %Finv = ; %rescaled and so that we exploit curvature of F_Xinv
        x = F_Xinv(P*(F_X(I(2)) -  F_X(I(1)))+ F_X(I(1)));
    end
end