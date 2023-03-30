function plot_from_coef_star(cn_center,cn,xq,yq,priorpar)

%N = 512;
%angles = linspace(0,2*pi,N)';
n = sqrt(length(xq));

%evalpar = makeEvalMatrixFourier(priorpar.angles,priorpar.maxfreq);
% evaluate Fourier basis in those points given coefficients cn
thetaq = priorpar.B * (cn.*priorpar.lambdahalf);
thetaq = priorpar.mean + thetaq;
% push-forward
gamma = push_forward_star2D(cn_center,thetaq,xq,yq,priorpar);
% reshape xq and yq coming from meshgrid
X = reshape(xq,n,n);
Y = reshape(yq,n,n);
Z = reshape(gamma,n,n);
% plot
surf(X,Y,Z,'EdgeColor','none','FaceColor','interp')
view(2)