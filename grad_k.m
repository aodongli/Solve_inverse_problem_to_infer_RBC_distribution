function [ grad ] = grad_k( k, dfdk_com, f, G_train, a1c_train, b )
%dfdk_com = @(t,k) (supp*exp(-k*supp)*(exp(-k*supp) - exp(-k*t)))./(exp(-k*supp) - 1)^2 - (supp*exp(-k*supp) - t.*exp(-k*t))./(exp(-k*supp) - 1);
dfdk = deriv_k(dfdk_com, b, k);

%f = @(k) (exp(-k*b) - exp(-k*supp)) / (1-exp(-k*supp));
grad = (G_train*f(k)-a1c_train)'*G_train*dfdk;


end

