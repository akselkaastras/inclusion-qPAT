function int = triangint3grad(g,kappa)
  
  
  % The function int=grinprodgaus(g,il,im); calculates the gradient part
  % in the linear  FEM in Optical Tomography
  % P. Ronkanen and M. Vauhkonen 10.5. 1996
  % Modified by V.Kolehmainen 9.7.1998
  %
  % Modified by T. Vilhunen 27.9.2001  
  %
  % Modified by N. Hänninen June 2017: 
  % Calculates the gradient part of FEM-matrix when kappa is in piecewise
  % linear basis.
  
    
w = [1/6*ones(3,1)];
ip = [1/2 0;1/2 1/2;0 1/2];
L = [-1 1 0;-1 0 1];
Jt = L*g;
iJt = inv(Jt);
dJt = abs(det(Jt));
G = iJt*L;
int = 0;
for ii = 1:3,
  S = [1-ip(ii,1)-ip(ii,2);ip(ii,1);ip(ii,2)];
  int_1 =  kappa(1)*S(1)*w(ii)*G'*G;
  int_2 =  kappa(2)*S(2)*w(ii)*G'*G;
  int_3 =  kappa(3)*S(3)*w(ii)*G'*G;
  int = int + int_1 + int_2 + int_3;
end
int = int*dJt;

