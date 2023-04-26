function plot_from_coef_2D(cn,priorpar)
eps = 0.1;

[X,Y] = meshgrid(linspace(-1-eps,1+eps,800));
xq = X(:);
yq = Y(:);
N = length(xq);
n = floor(sqrt(N));

% make evaluation matrix for Fourierbasis in query points (xq,yq)
evalpar = makeEvalMatrixFourier_2d(xq,yq,priorpar.maxfreq);
% evaluate Fourier basis in those points given coefficients cn
thetaq = evalpar.B * (cn.*priorpar.lambdahalf);
% reshape xq and yq coming from meshgrid
X = reshape(xq,n,n);
Y = reshape(yq,n,n);
Z = reshape(thetaq,n,n);
t = linspace(0,2*pi,10000);
x = cos(t);
y = sin(t);
z = 0*t+100;
%ind = X.^2 + Y.^2 > 1;
%Z(ind) = nan;
% plot
surf(X,Y,Z,'EdgeColor','none','FaceColor','interp')
hold on
plot3(x,y,z,'k-','linewidth',2);
view(2)
box off
axis off
axis equal
xlim([-1.1 1.1]) 
ylim([-1.1 1.1])