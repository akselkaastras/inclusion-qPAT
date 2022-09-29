function int = triangint3gradPre(g, quadPoints, k)


% The function int=grinprodgaus(g,il,im); calculates the gradient part
% in the linear  FEM in Optical Tomography
% P. Ronkanen and M. Vauhkonen 10.5. 1996
% Modified by V.Kolehmainen 9.7.1998
%
% Modified by T. Vilhunen 27.9.2001
%
% Modified by N. H�nninen June 2017:
% Calculates the gradient part of FEM-matrix in piecewise linear base.

switch quadPoints
    case 3
        w = 1/6*ones(3,1);
        ip = [1/2 0;1/2 1/2;0 1/2];
    case 7
        w = [1/40*ones(3,1); 1/15*ones(3,1);27/120];
        ip = [0 0; 1 0; 0 1; 1/2 0; 1/2 1/2; 0 1/2; 1/3 1/3];
end


w = [1/6*ones(3,1)];
ip = [1/2 0;1/2 1/2;0 1/2];
L = [-1 1 0;-1 0 1];
Jt = L*g;
iJt = inv(Jt);
dJt = abs(det(Jt));
G = iJt*L;
int = 0;
for ii = 1:length(w)
    S = [1-ip(ii,1)-ip(ii,2);ip(ii,1);ip(ii,2)];
    int = int + S(k)*w(ii)*G'*G;
end
int = int*dJt;

