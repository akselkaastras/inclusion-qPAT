function plot_from_coef_star_mean(results,priorpar)

xi_center = priorpar.center;
[X,Y] = meshgrid(linspace(-1,1,200));
xq = X(:);
yq = Y(:);
n = sqrt(length(xq));

M = size(results.XR_burnthin,1);
gamma = zeros(length(xq),1);
for i = 1:M
    xi = results.XR_burnthin(i,:);
    xi = reshape(xi,priorpar.dim);
    % evaluate Fourier basis in those points given coefficients cn
    thetaq = priorpar.B * xi;
    thetaq = priorpar.mean + thetaq;
    % push-forward
    gamma = gamma + push_forward_star2D(xi_center,thetaq,xq,yq,priorpar);
    i
end
gamma = gamma/M;
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