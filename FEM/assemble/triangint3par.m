function int = triangint3par(g, quadPoints, k, kappa)
  

switch quadPoints
case 3
  w = [1/6*ones(3,1)];
  ip = [1/2 0;1/2 1/2;0 1/2];
case 7
  w = [1/40*ones(3,1); 1/15*ones(3,1);27/120];
  ip = [0 0; 1 0; 0 1; 1/2 0; 1/2 1/2; 0 1/2; 1/3 1/3];
end
quad = g(1,:)+ip(:,1)*(g(2,:)-g(1,:))+ip(:,2)*(g(3,:)-g(1,:));
L = [-1 1 0;-1 0 1];
Jt = L*g;
dJt = abs(det(Jt));
kap = kappa(quad(:,1),quad(:,2));

int = 0;

for ii = 1:length(w)
  S = [1-ip(ii,1)-ip(ii,2);ip(ii,1);ip(ii,2)];
  int = int + w(ii)*S(k)*(kap(ii)*S');
end

int = int*dJt;
