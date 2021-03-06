function [muAS, vAS] = Antithetic(phi, Finv, N)

U = rand(N,1); %uniform distribustion

V = phi(Finv(U));
Vt = phi(Finv(1-U));

mu1 = (1/N)*sum(V);
mu2 = (1/N)*sum(Vt);

muAS = (mu1+mu2)/2;

vAS = (var(U)+cov(U,1-U))/2;

end