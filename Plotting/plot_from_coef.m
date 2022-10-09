function plot_from_coef(cn,xq,yq,priorpar)

N = length(xq);
n = floor(sqrt(N));

% make evaluation matrix for Fourierbasis in query points (xq,yq)
evalpar = makeEvalMatrixFourier(xq,yq,priorpar.maxfreq);
% evaluate Fourier basis in those points given coefficients cn
thetaq = evalKLFourierExpansion(cn,evalpar,priorpar);
% push-forward
gamma = push_forward_levelset2D(thetaq,priorpar);
% reshape xq and yq coming from meshgrid
X = reshape(xq,n,n);
Y = reshape(yq,n,n);
Z = reshape(gamma,n,n);
% plot
surf(X,Y,Z,'EdgeColor','none','FaceColor','interp')
view(2)