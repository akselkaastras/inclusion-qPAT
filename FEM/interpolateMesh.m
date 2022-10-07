function vq = interpolateMesh(u,xq,yq,meshpar)

%% Load mesh
p = meshpar.p';
H = meshpar.t(1:3,:)';

pN = size(p,1);
HN = size(H,1);

vq = zeros(length(xq),1);

z = 1:length(xq);
%% Run trough each element
for j = 1:HN
    ind = H(j,:);
    gg = p(ind,:);
    [IN, ON] = inpolygon(xq,yq,gg(:,1),gg(:,2));
    IN = bitor(IN,ON);
    if any(IN)
        id = z(IN);
        abc = [gg u(ind)];
        n = cross(abc(1,:)-abc(2,:),abc(3,:)-abc(2,:));
        n = n/n(3);
        vq(id) = -n(1)*(xq(id)-gg(2,1)) - n(2)*(yq(id)-gg(2,2)) + u(ind(2));
    else
    end
end