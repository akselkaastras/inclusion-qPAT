function int = triangint3gradpar(g,kappa)
  

    
w = [1/6*ones(3,1)];
ip = [1/2 0;1/2 1/2;0 1/2];
L = [-1 1 0;-1 0 1];
Jt = L*g;
iJt = inv(Jt);
dJt = abs(det(Jt));
G = iJt*L; % gradient of the three local basis functions in the "element space"
G = G';
int = 0;

S = [1-ip(1,1)-ip(1,2),1-ip(2,1)-ip(2,2),1-ip(3,1)-ip(3,2);ip(1,1),ip(2,1),ip(3,1);ip(1,2),ip(2,2),ip(3,2)];
%%
int = 0;
for ii = 1:3,

  int_1 =  kappa(1)*w(ii)*S(:,ii)*G(1,:)*G';
  int_2 =  kappa(2)*w(ii)*S(:,ii)*G(2,:)*G';
  int_3 =  kappa(3)*w(ii)*S(:,ii)*G(3,:)*G';
  int = int + int_1 + int_2 + int_3;
end
int = int*dJt;

