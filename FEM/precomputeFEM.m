function fmdl = precomputeFEM(meshpar,Ntrig)
%   Computes and saves integrals for FEM equations for cond eq.
%   Niko Hänninen 2019
%   Aksel Rasmussen 2021

% $$$ keyboard

p = meshpar.p';
H = meshpar.t(1:3,:)';
Ebound = meshpar.btri_ind;
Nbound = meshpar.e(1,:)';

pN = size(p,1);
HN = size(H,1);

Agrad = sparse(pN*pN,pN);

for kk = 1:pN
    if rem(kk,round(pN/20)) == 0
        disp(['- - ', num2str(round((kk/pN)*100)),' %'])
    end
    K = sparse(pN,pN);
    nodeInd = kk;
    for jj = 1:3
        elInd = find(H(:,jj) == nodeInd);
        for ii = elInd'
            ind = H(ii,:);
            gg = p(ind,:);
            elIntGrad = triangint3gradPre(gg, 3, jj);
            K(ind, ind) = K(ind, ind) + elIntGrad;
        end
    end
    Agrad(:,kk) = K(:);
end


Q = sparse(pN,2*Ntrig);
Nvec  = [[-Ntrig : -1],[1 : Ntrig]];

disp('- Building rhs')
for nn = 1:length(Nvec)
    n = Nvec(nn);
    for ii = 1:HN
        ind = H(ii,:);
        e1 = find(Ebound==ii,1);
        if ~isempty(e1)
            a = ismember(ind(:),Nbound);
            a1 = find(a==1);
            if length(a1)==2
                bind = ind(a1);
                bgg = p(bind,:);
                bint = boundint2phi(bgg,n);
                Q(bind,nn) = Q(bind,nn) + bint;
            else
                disp('something wrong')
                keyboard
            end
        end 

    end
end

% Insert inclusion in (0,0)
fcn = @(x,y) abs(x+1i*y)<0.3;
cond = 1+feval(fcn,meshpar.p(1,:),meshpar.p(2,:))';
N = length(cond);
condmatrix = spdiags(cond, 0, N, N);
K = Agrad*condmatrix;
K = reshape(sum(K,2),N,N);
K = 1/2*(K'+K);
K = K+1e-12*speye(size(K));
% Optimal ordering
[phi,stats] = amd(K); 
r(phi) = 1:max(size(phi));
Qperm = Q(phi,:);
fmdl.Agrad = Agrad;
fmdl.Q = Q;
fmdl.Qperm = Qperm;
fmdl.phi = phi;
fmdl.r = r;
fmdl.Nvec = Nvec;
fmdl.Ntrig = Ntrig;
fmdl.I = 1e-12*speye(size(K));

end

