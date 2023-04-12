%% Test precomputeRHS
meshpar_fine = mesh_comp(3);
sigma = 0.8;
m = 0.5*[sqrt(2),sqrt(2)];
scale = 4;

% Source is Gaussian
%wfun = @(x1,x2) scale*exp(-1/(2*sigma^2)*((x1-m(1)).^2+(x2-m(2)).^2));
%wfungrad = @(x1,x2) 1/(sigma^2)*wfun(x1,x2).*[m(1)-x1 m(2)-x2]; 
wfun = @(x1,x2) 2*x1.*x2+3;
wfungrad = @(x1,x2) [x1 x2];

% Precomputing finite element matrices and rhs
fmdl = precomputeFEM(meshpar_fine);
fmdl = precomputeRHStest(meshpar_fine,fmdl,wfun,wfungrad);
%% Testing Q1

pN = length(meshpar_fine.p);
D = sqrt(meshpar_fine.p(1,:).^2 + meshpar_fine.p(2,:).^2);
Dmatrix = spdiags(D', 0, pN, pN);

Q1 = fmdl.L1*Dmatrix;
Q1 = sum(Q1,2);
v = 0.5*meshpar_fine.p(1,:).^2 + 0.5*meshpar_fine.p(2,:).^2;
v*Q1

%% Testing Q2
D = meshpar_fine.p(1,:).^2 + meshpar_fine.p(2,:).^2;
Dmatrix = spdiags(D', 0, pN, pN);
Q2 = fmdl.L2*Dmatrix;
Q2 = sum(Q2,2);
v = 0*meshpar_fine.p(1,:)+1;
v*Q2
fun = @(r,theta) r.^2.*(2*r.*cos(theta).*r.*sin(theta)+3).*r;
q = integral2(fun,0,1,0,2*pi)