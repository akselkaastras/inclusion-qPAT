function plot_from_coef_star(xi,priorpar)

xi_center = priorpar.center;
[X,Y] = meshgrid(linspace(-1,1,200));
xq = X(:);
yq = Y(:);
n = sqrt(length(xq));

% evaluate Fourier basis in those points given coefficients cn
thetaq = priorpar.B * (xi.*priorpar.lambdahalf);
thetaq = priorpar.mean + thetaq;
% push-forward
gamma = push_forward_star2D(xi_center,thetaq,xq,yq,priorpar);
% reshape xq and yq coming from meshgrid
X = reshape(xq,n,n);
Y = reshape(yq,n,n);
Z = reshape(gamma,n,n);
% plot
figure(1);
surf(X,Y,Z,'EdgeColor','none','FaceColor','interp')
view(2)