%% Test precomputeFEM originally by Niko
clc; clear;
addpath(genpath(pwd))
%%
meshpar = mesh_comp(2);
mesh.r = meshpar.p';
mesh.H = meshpar.t(1:3,:)';
%%
rN = size(mesh.r,1);
HN = size(mesh.H,1);

Aint = sparse(rN*rN,rN);
Agrad = sparse(rN*rN,rN);
disp('- Building C and K')


for kk = 1:rN
    if rem(kk,round(rN/20)) == 0
        disp(['- - ', num2str(round((kk/rN)*100)),' %'])
    end
    C = sparse(rN,rN);
    K = sparse(rN,rN);
    nodeInd = kk;
    for jj = 1:3
        elInd = find(mesh.H(:,jj) == nodeInd);
        for ii = elInd'
            ind = mesh.H(ii,:);
            gg = mesh.r(ind,:);
            elInt = triangint3Pre(gg, 7, jj);
            elIntGrad = triangint3gradPre(gg, 3, jj);
            C(ind, ind) = C(ind, ind) + elInt;
            K(ind, ind) = K(ind, ind) + elIntGrad;
        end
    end
    Aint(:,kk) = C(:);
    Agrad(:,kk) = K(:);
    % Ct = Ct + C;
    % Kt = Kt + K;
end
%% What if I instead just use precomputeFEM out the bag
fmdl = precomputeFEM(meshpar);

%%
gamfun = @(x,y) 3*x.^2-y+2;
gamma = gamfun(mesh.r(:,1),mesh.r(:,2));
v = ones(rN,1);
gammatrix = spdiags(v, 0, rN, rN);

% Build mass matrix
C = Aint*gammatrix;
C = reshape(sum(C,2),rN,rN);
C = 1/2*(C+C');

% Build mass matrix for precomputeFEM
%C2 = fmdl.Aarea*gammatrix;
%C2 = reshape(sum(C2,2),rN,rN);
%C2 = 1/2*(C2+C2');
%% 
gamma'*C*gamma
gamma'*fmdl.Carea*gamma
fun = @(r,theta) r.*((3*abs(r.*cos(theta)).^2 - (r.*sin(theta))+2)).^2;
normsq = integral2(fun,0,1,0,2*pi)