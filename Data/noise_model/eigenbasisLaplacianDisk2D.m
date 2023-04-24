function [E,Lambda] = eigenbasisLaplacianDisk2D(MBessel,Nzeros,meshpar)

hmax = meshpar.hmax;



filename1 = ['eigenf_',num2str(hmax),'_',num2str(MBessel),'_',num2str(Nzeros),'_large.mat'];
filename2 = ['eigenv_',num2str(hmax),'_',num2str(MBessel),'_',num2str(Nzeros),'_large.mat'];

if ~isfile(['Data/noise_model/',filename1])
    % Mesh in polar coordinates
    [theta,r] = cart2pol(meshpar.p(1,:),meshpar.p(2,:));
    m = (0:MBessel)';
    E = zeros(2*MBessel+1,Nzeros,length(r));
    Lambda = zeros(2*MBessel+1,Nzeros);

    % Find N first Bessel zeros
    kmn = besselzero(m, Nzeros, 1);
    %m = (1:MBessel)';
    
    
    for mm = 1:MBessel
        x = kmn(mm+1,:)' * r;
        %J(mm+1,:,:) = besselj(mm,x);
        J = besselj(mm,x);
        nrmlz = sqrt(2)./(sqrt(pi)*abs(besselj(mm+1, kmn(mm+1,:))))';
        E(mm+1,:,:) = nrmlz.*(J.*cos(mm*theta));
        E(end-mm+1,:,:) = nrmlz.*(J.*sin(mm*theta));
        Lambda(mm+1,:) = kmn(mm+1,:);
        Lambda(end-mm+1,:) = kmn(mm+1,:);
        disp([' - - ',num2str(mm),'/',num2str(MBessel)])
    end
    % L = mm
    % Jzero = kmn(mm+1,n)
    

    nrmlz = sqrt(2)./(sqrt(2*pi)*abs(besselj(1, kmn(1,:))))';
    E(1,:,:) = nrmlz.*besselj(0,kmn(1,:)' * r);
    Lambda(1,:) = kmn(1,:);

    E = reshape(E,(2*MBessel+1)*Nzeros,length(r));
    Lambda = reshape(Lambda,(2*MBessel+1)*Nzeros,1);
    Lambda = Lambda.^2;

    save(['Data/noise_model/',filename1],'E','-v7.3')
    save(['Data/noise_model/',filename2],'Lambda','-v7.3')
else
    s1 = load(['Data/noise_model/',filename1]);
    s2 = load(['Data/noise_model/',filename2]);
    E = s1.E;
    Lambda = s2.Lambda;
    disp('Eigenfunctions and eigenvalues loaded')
end

