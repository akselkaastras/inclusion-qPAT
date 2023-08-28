function vq = interpolateMesh(u,xq,yq,meshpar)
% Returns vq as u evaluated in query points (xq,yq) based on the function u
% represented in a mesh specified by meshpar

%% Load mesh
p = meshpar.p';
H = meshpar.t(1:3,:)';

pN = size(p,1);
HN = size(H,1);

vq = zeros(length(xq),1);

z = 1:length(xq);


%% Are there points in (xq,yq) outside the mesh meshpar_fine?
% Only works if domain is convex
xv = meshpar.p(1,:)';
yv = meshpar.p(2,:)';
kk = convhull(xv,yv);
%[in, on] = inpolygon(xq,yq,xv(k),yv(k));
[in,on] = inpoly2([xq yq],[xv(kk) yv(kk)],[],1e-7);
indzero = ~(in+on);

ind = 1:pN;
p0 = meshpar.p(:,meshpar.e(1,:));
p1 = meshpar.p(:,meshpar.e(2,:));

% project to nearest line
for k = ind(indzero)
    dist = distToLine(p0,p1,[xq(k),yq(k)]);
    [~,I] = min(dist);
    znew = proj(p0(:,I), p1(:,I), [xq(k),yq(k)]);
    xq(k) = znew(1);
    yq(k) = znew(2);
end

[in,on] = inpoly2([xq yq],[xv(kk) yv(kk)],[],1e-7);
indzero = ~(in+on);
%disp([num2str(sum(indzero)),' points outside mesh']);

%% Run through each element
for j = 1:HN
    ind = H(j,:);
    gg = p(ind,:);
    [IN, ON] = inpoly2([xq, yq],[gg(:,1), gg(:,2)],[],1e-4);
    IN = bitor(IN,ON);
    if any(IN)
        id = z(IN);
        abc = [gg u(ind)];
        n = cross(abc(1,:)-abc(2,:),abc(3,:)-abc(2,:));
        n = n/n(3);
        vq(id) = -n(1)*(xq(id)-gg(2,1)) - n(2)*(yq(id)-gg(2,2)) + u(ind(2));
        if abs(vq(id)) < 1e-8
            %disp('small value');
            %keyboard;
        end
    else
    end
end

