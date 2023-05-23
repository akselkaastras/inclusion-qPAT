function [A,B,C] = computeProjMatrix(meshpar)

%% Load mesh
p = meshpar.p';
H = meshpar.t(1:3,:)';
Ebound = meshpar.btri_ind;
Nbound = meshpar.e(1,:)';

pN = size(p,1);
NN = size(Nbound,1);
HN = size(H,1);

pNN = pN-NN;
%% Matrices
A = sparse(HN,pN);
B = sparse(HN,pN);
C = sparse(HN,pN);
kappaA = @(x,y) x;
kappaB = @(x,y) y;
kappaC = @(x,y) 0*x+1;

for ii = 1:HN
    if rem(ii,round(HN/20)) == 0
        disp(['- - Assembling A, B and C - - ', num2str(round((ii/HN)*100)),' %'])
    end
    ind = H(ii,:);
    gg = p(ind,:);
    elIntA = triangint3par1(gg, 3, kappaA);
    elIntB = triangint3par1(gg, 3, kappaB);
    elIntC = triangint3par1(gg, 3, kappaC);
    A(ii, ind) = A(ii, ind) + elIntA;
    B(ii, ind) = B(ii, ind) + elIntB;
    C(ii, ind) = C(ii, ind) + elIntC;
end
