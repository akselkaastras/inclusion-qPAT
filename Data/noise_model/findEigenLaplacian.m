function findEigenLaplacian(meshpar)

% 4 GB for each core
%LL = 10;
LL = 500;
[K,C,M] = computeLaplacian(meshpar);
%E = zeros(size(K,1),L);
%D = zeros(L,1);
NZ = setdiff(1:length(meshpar.p),meshpar.e(1,:));
pN = length(meshpar.p);

filename = 'Data/noise_model/eigenv/eigen_0.01';

parfor i = 1:30
    
    [VV,D,iresult] = sptarn(K,C,LL*(i-1),LL*i,1,1e-13,1000);
    V = zeros(pN,size(VV,2));
    for j = 1:size(VV,2)
        V(NZ,j) = VV(:,j);
        V(:,j) = V(:,j)/sqrt(V(:,j)'*M*V(:,j));
    end
    parsave(strcat(filename,'_',num2str(i),'.mat'),V,D);
    if iresult <0
        disp('Did not find all eigenvalues')
    end

end