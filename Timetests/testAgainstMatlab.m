clc;
clear;
%% Solve pde in MATLAB
model = createpde();
geometryFromEdges(model,@circleg);
hmax = 0.01;
generateMesh(model,"Hmax",hmax,'GeometricOrder','linear');
pdemesh(model); 
meshpar = mesh_comp(hmax);

% define source from wfun
sigma = 0.5;
m = 0.5*[sqrt(2),sqrt(2)];
scale = 10;
figure(2);

myFun = @(location,state) scale*exp(-1/(2*sigma^2)*((location.x-m(1)).^2+(location.y-m(2)).^2)); 
% Source is Gaussian
wfun = @(x1,x2) scale*exp(-1/(2*sigma^2)*((x1-m(1)).^2+(x2-m(2)).^2));
wfungrad = @(x1,x2) 1/(sigma^2)*wfun(x1,x2).*[m(1)-x1 m(2)-x2]; 

a = @(location,state)100*location.y.^2.*tanh(location.x).^2;
c = @(location,state) 2+10*((location.x+0.3).^2 + location.y.^2).^(1/2);

applyBoundaryCondition(model,"dirichlet", ...
                             "Edge",1:model.Geometry.NumEdges, ...
                             "u",myFun);

specifyCoefficients(model,"m",0,"d",0,"c",c,"a",a,"f",0);



results = solvepde(model);

u_matlab = results.NodalSolution;
pdeplot(model,"XYData",u_matlab)
title("Numerical Solution");
xlabel("x")
ylabel("y")

%% Solve PDE with own code
location.x = meshpar.p(1,:);
location.y = meshpar.p(2,:);
D = c(location,0);
gamma = a(location,0);
% Precomputing finite element matrices and rhs
fmdl = precomputeFEM(meshpar);
fmdl = precomputeRHS(meshpar,fmdl,wfun,wfungrad);
% Precomputing stiffness

fmdl = fixingD(meshpar,fmdl,D);
% Evaluating forward model from precomputed matrices
u = evalFowardModel(fmdl,meshpar,gamma');
% Add source back again
pN = length(meshpar.p);
meshpar.NZ = setdiff(1:length(meshpar.p),meshpar.e(1,:));

U = zeros(pN,1);
U(meshpar.NZ) = u;

u = U + wfun(meshpar.p(1,:)',meshpar.p(2,:)');
figure;
plot_from_gamma(U,meshpar)
norm(u_matlab-u,2)