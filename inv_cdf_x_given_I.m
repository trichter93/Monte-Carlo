function Finv = inv_cdf_x_given_I(U, F_X, I)
    if(U==0)
        Finv = I(1);
    elseif (U==1)
        Finv = I(2);
    else
        points=1e4;
        x = linspace(I(1), I(2), points);
        u = linspace(0, 1, points);
        dp=1e-5;
        for i=1:points
            if(abs(U-cdf_x_given_I(x(i), F_X, I))<=dp)
                Finv = (u(i+1) - U)*x(i)/(u(i+1)-u(i)) + (U - u(i))*x(i+1)/(u(i+1)-u(i));
                return;
            end
        end
    end
end