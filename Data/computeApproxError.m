function datapar = computeApproxError(datapar,N,priortype)

% Input:
%   datapar: struct
%   N: number of iterations
%   priortype: priorpar.type ('level','star','starDG')
%
% Output:
%   datapar: updated datapar containing:
%       sigmasq: sigmasq>0 such that N(0,sigmasq I) is the best approximation 
%       of the empirical distribution of G_{fine}-G_{coarse}

%% Is the approximation error already computed?
filename = strcat(num2str(datapar.meshpar_fine.hmax),'_',num2str(datapar.meshpar.hmax),'_',priortype,'_',num2str(datapar.noiselevel),'_',num2str(datapar.N), '.mat');
disp(filename);
if ~isfile(strcat('Data/sigma2/sigma2_',filename))

    %% Initialize mesh
    
    meshpar_fine = datapar.meshpar_fine;
    meshpar = datapar.meshpar;
    pN_coarse = length(meshpar.p);
    pN = length(meshpar_fine.p);
    meshpar_fine.NZ = setdiff(1:pN,meshpar_fine.e(1,:));
    meshpar.NZ = setdiff(1:pN_coarse,meshpar.e(1,:));
    
    %% Initialize prior

    % Here we choose some prior parameters suitable for samples
    if strcmpi(priortype, 'level')
        q = 4;
        alpha = 1.2;
        tau = 1;
        maxfreq = 4;
        delta = 0.1;
        priorpar = prior_init_2d(meshpar_fine.p(1, :), meshpar_fine.p(2, :),alpha,tau,q,maxfreq);
        priorpar_coarse = prior_init_2d(meshpar.p(1, :), meshpar.p(2, :),alpha,tau,q,maxfreq);
    
        %priorpar = prior_init_2d(meshpar_fine.xq', meshpar_fine.yq',alpha,tau,q,maxfreq);
        %priorpar_coarse = prior_init_2d(meshpar.xq', meshpar.yq',alpha,tau,q,maxfreq);
    
        % For levelset
        priorpar.ninterface = 3; % number of interfaces
        priorpar.c = [-inf 0 2 inf]; % contour levels (ninterface+1)
        priorpar.v = [0.1 0.3 0.5]; % values at each interface
        priorpar.delta = delta; % smoothing factor
        priorpar.dim = [priorpar.M,1];
        priorpar.type = 'level';
        priorpar.std = 1;
        priorpar_coarse.ninterface = 3;
        priorpar_coarse.c = [-inf 0 2 inf]; % contour levels (ninterface+1)
        priorpar_coarse.v = [0.1 0.3 0.5]; % values at each interface
        priorpar_coarse.delta = delta; % smoothing factor
        priorpar_coarse.dim = [priorpar.M,1];
        priorpar_coarse.type = 'level';
        priorpar_coarse.std = 1;
    elseif strcmpi(priortype,'star')
        q = 10;
        alpha = 1;
        tau = 0.5;
        maxfreq = 7;
        xq = linspace(0,2*pi,256);
        
        priorpar = prior_init(xq,alpha,tau,q,maxfreq);
        
        % For star-shaped inclusion parametrization
        priorpar.ninclusions = 2; % number of interfaces
        priorpar.v = [0.2 0.4]; % values at each interface
        
        priorpar.background = 0.1;
        priorpar.angles = linspace(0,2*pi,256)';
        priorpar.mean = -2;
        priorpar.dim = [priorpar.M,priorpar.ninclusions];
        priorpar.type = 'star';
        priorpar.center = [0.37,-0.43;-0.44,0.36];
        priorpar.std = 0.2;
    elseif strcmpi(priortype,'starDG')
        q = 30;
        alpha = 1.75;
        tau = 1;
        maxfreq = 12;
        xq = linspace(0,2*pi,512);
        
        priorpar = prior_init(xq,alpha,tau,q,maxfreq);
        
        % For star-shaped inclusion parametrization
        priorpar.ninclusions = 2; % number of interfaces
        priorpar.v = [0.2 0.4]; % values at each interface
        
        priorpar.background = 0.1;
        priorpar.angles = linspace(0,2*pi,512)';
        priorpar.mean = -2;
        priorpar.dim = [priorpar.M,priorpar.ninclusions];
        priorpar.type = 'starDG';
        priorpar.center = [0.37,-0.43;-0.44,0.36];
        priorpar.std = 0.2;
    end

    %% Initialize forward model on fine mesh
    D = datapar.D;

    if strcmpi(priortype,'starDG')
        fmdl = precomputeFEM_DG(meshpar_fine);
        fmdl = precomputeRHS_DG(meshpar_fine,fmdl,datapar.wfun,datapar.wfungrad);
    else
        fmdl = precomputeFEM(meshpar_fine);
        fmdl = precomputeRHS(meshpar_fine,fmdl,datapar.wfun,datapar.wfungrad);
    end


    % Precomputing stiffness
    fmdl = fixingD(meshpar_fine,fmdl,D);

    M = datapar.N;
    trunc = M*(2*M+1);
    % Projection matrices for projection
    fmdl = computeProjectionMatrices_fine(fmdl,meshpar_fine,trunc);    
    %% Initialize forward model on coarse mesh
    
    D_coarse = interpolateMesh(D',meshpar.p(1,:)',meshpar.p(2,:)',meshpar_fine);
    

    %if strcmpi(priortype,'starDG') || strcmpi(priortype,'level')
    if strcmpi(priortype,'starDG')
        fmdl_coarse = precomputeFEM_DG(meshpar);
        fmdl_coarse = precomputeRHS_DG(meshpar,fmdl_coarse,datapar.wfun,datapar.wfungrad);
        % Precomputing stiffness
        fmdl_coarse = fixingD(meshpar,fmdl_coarse,D_coarse');

    else
        % Precomputing finite element matrices and rhs
        fmdl_coarse = precomputeFEM(meshpar);
        fmdl_coarse = precomputeRHS(meshpar,fmdl_coarse,datapar.wfun,datapar.wfungrad);
        
        % Precomputing stiffness
        fmdl_coarse = fixingD(meshpar,fmdl_coarse,D_coarse');
    end
   

    % Projection matrices for coarse mesh
    E = eigenbasisFEM(meshpar_fine,trunc);
    fmdl.U_proj = fmdl.Carea*E;
    fmdl_coarse = computeProjectionMatrices_coarse(fmdl_coarse,meshpar,datapar.meshpar_fine,priorpar,trunc);
    
    %% Compute samples of G_{fine}(prior_sample) - G_{coarse}(prior_sample)
    V = zeros(trunc,N);
    U = zeros(pN,1);
    disp('Computing samples of G_{fine}(prior_sample) - G_{coarse}(prior_sample) for approximation error')
    for i = 1:N
        disp(['Sample ',num2str(i),'/',num2str(N)])
        % Sample prior
        
        xi = priorpar.std*randn(priorpar.dim);
        xi = priorpar.lambdahalf.*xi;
        theta = priorsample(xi,priorpar);
        
        % Push-forward
        if strcmpi(priortype,'level')
            gamma = push_forward_levelset2D_smooth(theta,priorpar);
            theta_coarse = priorsample(xi,priorpar_coarse);
            gamma_coarse = push_forward_levelset2D_smooth(theta_coarse,priorpar);
            
            %gamma = push_forward_level2D_interp(theta,meshpar_fine.n,meshpar_fine.HN,priorpar);
            %theta_coarse = priorsample(xi,priorpar_coarse);
            %gamma_coarse = push_forward_level2D_interp(theta_coarse,meshpar.n,meshpar.HN,priorpar_coarse);
        elseif strcmpi(priortype,'star')
            xi_center = priorpar.center;
            theta = priorpar.mean+theta;
            % push-forward
            gamma = push_forward_star2D(xi_center,theta,datapar.meshpar_fine.p(1,:)',datapar.meshpar_fine.p(2,:)',priorpar);
            gamma_coarse = push_forward_star2D(xi_center,theta,datapar.meshpar.p(1,:)',datapar.meshpar.p(2,:)',priorpar);
        elseif strcmpi(priortype,'starDG')
            xi_center = priorpar.center;
            theta = priorpar.mean+theta;
            % push-forward
            gamma = push_forward_star2D_interp(xi_center,theta,meshpar_fine.xq,meshpar_fine.yq,meshpar_fine.n,meshpar_fine.HN,priorpar);
            gamma_coarse = push_forward_star2D_interp(xi_center,theta,datapar.meshpar.xq,datapar.meshpar.yq,datapar.meshpar.n,datapar.meshpar.HN,priorpar);
        end
        
        % Evaluating forward model from precomputed matrices
        u = evalFowardModel(fmdl,meshpar_fine,gamma);
        U(meshpar_fine.NZ) = u;
        u = U + datapar.wfun(meshpar_fine.p(1,:)',meshpar_fine.p(2,:)');
    

        u_coarse = evalFowardModel(fmdl_coarse,meshpar,gamma_coarse);
        U_coarse = zeros(length(meshpar.p),1);
        U_coarse(meshpar.NZ) = u_coarse;
        u_coarse = U_coarse + datapar.wfun(meshpar.p(1,:)',meshpar.p(2,:)');
        
        
        %if strcmpi(priortype,'starDG') || strcmpi(priortype,'level')
        if strcmpi(priortype,'starDG')
            % fine mesh
            u = u(meshpar_fine.t(1:3,:));
            u = fmdl.G*u(:);
            u = reshape(u,[3 meshpar_fine.HN]);

            data = gamma'.*u;
            data_proj = data(1,:)*fmdl.U_proj1 + data(2,:)*fmdl.U_proj2 + data(3,:)*fmdl.U_proj3;
            
            
            
            % coarse mesh
            u_coarse = u_coarse(meshpar.t(1:3,:));
            u_coarse = fmdl_coarse.G*u_coarse(:);
            u_coarse = reshape(u_coarse,[3 meshpar.HN]);
        
        
            data_coarse = gamma_coarse'.*u_coarse;
            data_proj_coarse = data_coarse(1,:)*fmdl_coarse.U_proj_coarse1 + data_coarse(2,:)*fmdl_coarse.U_proj_coarse2 + data_coarse(3,:)*fmdl_coarse.U_proj_coarse3;
       else
            % Multiply with gamma
            data = u.*gamma;
            data_coarse = u_coarse.*gamma_coarse;

            % Project to subspace
            data_proj = data'*fmdl.U_proj;
            data_proj_coarse = data_coarse'*fmdl_coarse.U_proj_coarse;
        end
        

        V(:,i) = data_proj'-data_proj_coarse';
    end
    
    %% Sample statistics
    m = mean(V,2);
    C = cov(V');
    sigmasq = (norm(m,2)^2+trace(C))/(length(m));
    eigC = eig(C);
    indpos = find(eigC>0,1);
    sigmasq_ws = (sum(eigC(indpos:end).^(1/2))/length(m))^2;
    sigmaKL2 = length(m)/sum(eigC(indpos:end).^(-1));
    sigmasq_m = trace(C)/length(m);
    sigma.sigmasq = sigmasq;
    sigma.sigmasq_ws = sigmasq_ws;
    sigma.sigmaKL2 = sigmaKL2;
    sigma.V = V;
    sigma.m = m;
    sigma.sigmasq_m = sigmasq_m;
    save(strcat('Data/sigma2/sigma2_',filename),'sigma')
else
    s = load(strcat('Data/sigma2/sigma2_',filename));
    sigmasq = s.sigma.sigmasq;
    sigmasq_ws = s.sigma.sigmasq_ws;
    sigmaKL2 = s.sigma.sigmaKL2;
    V = s.sigma.V;
    m = s.sigma.m;
    sigmasq_m = s.sigma.sigmasq_m;
    disp('Approximation error loaded')
end
datapar.filename = filename;
datapar.sigmasq = sigmasq;
datapar.sigmasq_ws = sigmasq_ws;
datapar.sigmaKL2 = sigmaKL2;
datapar.m = m;
datapar.sigmasq_m = sigmasq_m;
datapar.epssq_approx = datapar.epssq + datapar.sigmasq;
datapar.V = V;

