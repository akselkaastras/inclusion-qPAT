function error = computeL2error_level(results,priorpar,res)

[X,Y] = meshgrid(linspace(-1,1,res));
xq = X(:);
yq = Y(:);
n = sqrt(length(xq));

M = size(results.XR_burnthin,1);
gamma = zeros(length(xq),1);
% make evaluation matrix for Fourierbasis in query points (xq,yq)
evalpar = makeEvalMatrixFourier_2d(xq,yq,priorpar.maxfreq);

for i = 1:M
    xi = results.XR_burnthin(i,:);
    % evaluate Fourier basis in those points given coefficients cn
    thetaq = evalpar.B * xi';
    % push-forward
    gamma = gamma + push_forward_levelset2D_smooth(thetaq,priorpar);
    i
end
gamma = gamma/M;
% reshape xq and yq coming from meshgrid
X = reshape(xq,n,n);
Y = reshape(yq,n,n);
Z = reshape(gamma,n,n);

%% Calculate L^2 error 
% True phantom is defined as
% Make curves
t = linspace(0,2*pi,1000);
kitecurve = [cos(t)+0.65*cos(2*t)-0.65; 1.5*sin(t)];
r = sqrt(0.8+0.8*(cos(4*t)-1).^2);
cushioncurve = [r.*cos(t); r.*sin(t)]; 
curves = [0.18*kitecurve+[0.4;-0.4]; 0.12*cushioncurve+[-0.4;0.4]];
values = [0.2,0.4,0.1];

ninclusions = size(curves,1)/2;
ncurve = size(curves,2);

% Find mesh on boundary
npoints = length(xq);


% Make gamma from curves
HN = res*res;
Gamma = zeros(ninclusions,HN);
for i = 1:ninclusions
    index = 2*(i-1);
    nodes = [curves(index+1,:)',curves(index+2,:)'];
    Gamma(i,:) = Gamma(i,:) + inpoly2([xq, yq],nodes)';
end
gamma0 = values(ninclusions+1)+zeros(1,HN);
for i = 1:ninclusions
    gamma0 = gamma0 + values(i)*Gamma(i,:);
end

Z0 = reshape(gamma0,res,res);
ind = X.^2 + Y.^2 > 1;
Z(ind) = values(3);
area = 4/(res*res);
error = sqrt(area*sum(sum((Z-Z0).^2)));