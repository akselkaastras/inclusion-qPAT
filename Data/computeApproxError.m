function datapar = computeApproxError(datapar,N)

% Input:
%   datapar: struct
%   N: number of iterations
%
% Output:
%   sigma2: sigma2>0 such that N(0,sigma2 I) is the best approximation of 
%   the distribution of G_{fine}-G_{coarse}

%% Is the approximation error already computed?
filename = ['sigma2_',num2str(datapar.meshpar_fine.hmax),'_',num2str(datapar.meshpar.hmax),'.mat'];

if ~isfile(['Data/sigma2/',filename])

    %% Initialize mesh
    
    meshpar_fine = datapar.meshpar_fine;
    meshpar = datapar.meshpar;
    pN_coarse = length(meshpar.p);
    pN = length(meshpar_fine.p);
    meshpar_fine.NZ = setdiff(1:pN,meshpar_fine.e(1,:));
    meshpar.NZ = setdiff(1:pN_coarse,meshpar.e(1,:));
    
    %% Initialize forward model on fine mesh
    D = datapar.D;
    % Precomputing finite element matrices and rhs
    fmdl = precomputeFEM(meshpar_fine);
    fmdl = precomputeRHS(meshpar_fine,fmdl,datapar.wfun,datapar.wfungrad);
    
    % Precomputing stiffness
    fmdl = fixingD(meshpar_fine,fmdl,D);
    
    %% Initialize forward model on coarse mesh
    
    D_coarse = interpolateMesh(D',meshpar.p(1,:)',meshpar.p(2,:)',meshpar_fine);
    
    % Precomputing finite element matrices and rhs
    fmdl_coarse = precomputeFEM(meshpar);
    fmdl_coarse = precomputeRHS(meshpar,fmdl_coarse,datapar.wfun,datapar.wfungrad);
    
    % Precomputing stiffness
    fmdl_coarse = fixingD(meshpar,fmdl_coarse,D_coarse');
    
    %% Initialize prior
    q = 3;
    alpha = 1;
    tau = 1;
    maxfreq = 4;
    delta = 0.2;
    
    priorpar = prior_init_2d(meshpar_fine.p(1, :), meshpar_fine.p(2, :),alpha,tau,q,maxfreq);
    priorpar_coarse = prior_init_2d(meshpar.p(1, :), meshpar.p(2, :),alpha,tau,q,maxfreq);

    % For levelset
    priorpar.ninterface = 3; % number of interfaces
    priorpar.c = [-inf -2 2 inf]; % contour levels (ninterface+1)
    priorpar.v = [0.3 0.1 0.5]; % values at each interface
    priorpar.delta = delta; % smoothing factor
    priorpar.dim = [priorpar.M,1];
    priorpar.type = 'level';
    priorpar.std = 1;
    priorpar_coarse.ninterface = 3;
    priorpar_coarse.c = [-inf -2 2 inf]; % contour levels (ninterface+1)
    priorpar_coarse.v = [0.3 0.1 0.5]; % values at each interface
    priorpar_coarse.delta = delta; % smoothing factor
    priorpar_coarse.dim = [priorpar.M,1];
    priorpar_coarse.type = 'level';
    priorpar_coarse.std = 1;

    
    %% Compute samples of G_{fine}(prior_sample) - G_{coarse}(prior_sample)
    V = zeros(length(meshpar.p),N);
    U = zeros(pN,1);
    
    for i = 1:N
        i
        % Sample prior
        xi = randn(priorpar.dim);
        theta = priorsample(xi,priorpar);
        
        % Push-forward
        gamma = push_forward_levelset2D_smooth(theta,priorpar);
        gamma_coarse = push_forward_levelset2D_smooth(theta,priorpar_coarse);
        
        % Evaluating forward model from precomputed matrices
        u = evalFowardModel(fmdl,meshpar_fine,gamma);
        U(meshpar_fine.NZ) = u;
        u = U + datapar.wfun(meshpar_fine.p(1,:)',meshpar_fine.p(2,:)');
    
        uq = interpolateMesh(u,meshpar.p(1,:)',meshpar.p(2,:)',meshpar_fine);
        u_coarse = evalFowardModel(fmdl_coarse,meshpar,gamma_coarse);
        U_coarse = zeros(length(meshpar.p),1);
        U_coarse(meshpar.NZ) = u_coarse;
        u_coarse = U_coarse + datapar.wfun(meshpar.p(1,:)',meshpar.p(2,:)');
        %V(:,i) = uq(meshpar.NZ)-u_coarse;
        V(:,i) = uq-u_coarse;
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

