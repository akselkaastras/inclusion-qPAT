function int = triangint3gradpar(g, k, gradkappa)


w = [1/6*ones(3,1)];
ip = [1/2 0;1/2 1/2;0 1/2];
quad = g(1,:)+ip(:,1)*(g(2,:)-g(1,:))+ip(:,2)*(g(3,:)-g(1,:));

L = [-1 1 0;-1 0 1];
Jt = L*g;
iJt = inv(Jt);
dJt = abs(det(Jt));
G = iJt*L;
int = 0;
nabkap = gradkappa(quad(:,1),quad(:,2));

for ii = 1:length(w)
    S = [1-ip(ii,1)-ip(ii,2);ip(ii,1);ip(ii,2)];
    int = int + S(k)*w(ii)*nabkap(ii,:)*G;
end
int = int*dJt;