function Finv = inv_cdf_x_given_I(p, F_X, I)
    if(p==0)
        Finv = I(1);
    elseif (p==1)
        Finv = I(2);
    else
        points=1e4;
        x = linspace(I(1), I(2), points);
        dp=1e-5;
        for i=1:points
            if(p-cdf_x_given_I(x(i), F_X, I)<=dp)
                Finv = x(i);
                return;
            end
        end
    end
end