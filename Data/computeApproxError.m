function datapar = computeApproxError(datapar,N)

% Input:
%   datapar: struct
%   N: number of iterations
%
% Output:
%   sigma2: sigma2>0 such that N(0,sigma2 I) is the best approximation of 
%   the distribution of G_{fine}-G_{coarse}

%% Is the approximation error already computed?
filename = ['sigma2_',num2str(datapar.meshpar_fine.refine),'_',num2str(datapar.meshpar.refine),'.mat'];

if ~isfile(['Data/sigma2/',filename])

    %% Initialize mesh
    
    meshpar_fine = datapar.meshpar_fine;
    meshpar = datapar.meshpar;
    pN_coarse = length(meshpar.p);
    pN = length(meshpar_fine.p);
    meshpar_fine.NZ = setdiff(1:pN,meshpar_fine.e(1,:));
    meshpar.NZ = setdiff(1:pN_coarse,meshpar.e(1,:));
    
    %% Source function
    
    sigma = 0.5;
    m = 0.5*[sqrt(2),sqrt(2)];
    wfun = @(x1,x2) 2*exp(-1/(2*sigma^2)*((x1-m(1)).^2+(x2-m(2)).^2));
    wfungrad = @(x1,x2) 1/(sigma^2)*wfun(x1,x2).*[m(1)-x1 m(2)-x2]; 
    
    %% Initialize forward model on fine mesh
    D = datapar.D;
    % Precomputing finite element matrices and rhs
    fmdl = precomputeFEM(meshpar_fine);
    fmdl = precomputeRHS(meshpar_fine,fmdl,wfun,wfungrad);
    
    % Precomputing stiffness
    fmdl = fixingD(meshpar_fine,fmdl,D);
    
    %% Initialize forward model on coarse mesh
    
    D_coarse = interpolateMesh(D',meshpar.p(1,:)',meshpar.p(2,:)',meshpar_fine);
    
    % Precomputing finite element matrices and rhs
    fmdl_coarse = precomputeFEM(meshpar);
    fmdl_coarse = precomputeRHS(meshpar,fmdl_coarse,wfun,wfungrad);
    
    % Precomputing stiffness
    fmdl_coarse = fixingD(meshpar,fmdl_coarse,D_coarse');
    
    %% Initialize prior
    q = 3;
    alpha = 1;
    tau = 1;
    maxfreq = 3;
    delta = 0.2;
    
    priorpar = prior_init_2d(meshpar_fine.p(1, :), meshpar_fine.p(2, :),alpha,tau,q,maxfreq);
    priorpar_coarse = prior_init_2d(meshpar.p(1, :), meshpar.p(2, :),alpha,tau,q,maxfreq);
    
    % For levelset
    priorpar.ninterface = 3; % number of interfaces
    priorpar.c = [-inf -2 2 inf]; % contour levels (ninterface+1)
    priorpar.v = [1 0 2]; % values at each interface
    priorpar.delta = delta; % smoothing factor
    priorpar_coarse.ninterface = 3; % number of interfaces
    priorpar_coarse.c = [-inf -2 2 inf]; % contour levels (ninterface+1)
    priorpar_coarse.v = [1 0 2]; % values at each interface
    priorpar_coarse.delta = delta; % smoothing factor
    
    %% Compute samples of G_{fine}(prior_sample) - G_{coarse}(prior_sample)
    V = zeros(length(meshpar.NZ),N);
    U = zeros(pN,1);
    
    for i = 1:N
        i
        % Sample prior
        xi = randn(priorpar.M,1);
        theta = priorsample(xi,priorpar);
        
        % Push-forward
        gamma = push_forward_levelset2D_smooth(theta,priorpar);
        gamma_coarse = push_forward_levelset2D_smooth(theta,priorpar_coarse);
        
        % Evaluating forward model from precomputed matrices
        u = evalFowardModel(fmdl,meshpar_fine,gamma);
        U(meshpar_fine.NZ) = u;
        %u = U + wfun(meshpar.p(1,:)',meshpar.p(2,:)');
    
        uq = interpolateMesh(U,meshpar.p(1,:)',meshpar.p(2,:)',meshpar_fine);
        u_coarse = evalFowardModel(fmdl_coarse,meshpar,gamma_coarse);
        V(:,i) = uq(meshpar.NZ)-u_coarse;
    end
    
    %% Sample statistics
    m = mean(V,2);
    C = cov(V');
    sigmasq = (norm(m,2)^2+trace(C))/(length(m));
    save(['Data/sigma2/',filename],'sigmasq')
else
    s = load(['Data/sigma2/',filename]);
    sigmasq = s.sigmasq;
    disp('Approximation error loaded')
end
datapar.sigmasq = sigmasq;
datapar.epssq_approx = datapar.epssq + datapar.sigmasq;
