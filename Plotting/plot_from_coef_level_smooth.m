function plot_from_coef_level_smooth(cn,priorpar)


[X,Y] = meshgrid(linspace(-1,1,200));
xq = X(:);
yq = Y(:);
N = length(xq);
n = floor(sqrt(N));

% make evaluation matrix for Fourierbasis in query points (xq,yq)
evalpar = makeEvalMatrixFourier_2d(xq,yq,priorpar.maxfreq);
% evaluate Fourier basis in those points given coefficients cn
thetaq = evalpar.B * cn;
% push-forward
gamma = push_forward_levelset2D_smooth(thetaq,priorpar);
% reshape xq and yq coming from meshgrid
X = reshape(xq,n,n);
Y = reshape(yq,n,n);
Z = reshape(gamma,n,n);
ind = X.^2 + Y.^2 > 1;
Z(ind) = nan;
% plot
surf(X,Y,Z,'EdgeColor','none','FaceColor','interp')
view(2)
box off
axis off
axis equal
xlim([-1.1 1.1]) 
ylim([-1.1 1.1])