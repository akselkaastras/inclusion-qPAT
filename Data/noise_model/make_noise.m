function [xi, norm_noise, norm_noise_fem, wk] = make_noise(meshpar,fmdl,E,Lambda,trunc)

pN = length(meshpar.p);
if pN <= length(Lambda)
else
    error('Please choose more basisfunctions');
end
[Lambda,ord] = sort(Lambda);
E = E(ord,:);
E = E(1:trunc,:);
Lambda = Lambda(1:trunc);

e = randn(1,size(E,1));
xi = e * E;

% scale coefficients
wk = Lambda.^(-1);
scale_e = e'.*(wk);

% We could compute a different norm \| . \|_{-}, but we use L^2 for now 
%xi_smooth = scale_e' * E;
%plot_from_gamma(xi,meshpar);

norm_noise = sum(e.^2);
norm_noise_fem = xi*fmdl.Carea*xi';
