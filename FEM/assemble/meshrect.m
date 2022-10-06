function mesh = meshrect(a,b,c,d,M)
% returns struct with uniform mesh on domain
%% Generate mesh
%meshName = '15x10A0';
%mesh = load(['./Meshes/', meshName]);

% mesh.p % 36 points (36 x 2)
% mesh.H % element table (mesh: 1 = point 2,1,3) EtoV
% mesh.BH % line segments (node to node) that are in the boundary
% mesh.Ebound % which elements are at the boundary
% mesh.Domain % 1 if element is in domain
% mesh.Bdomain % 1 if element is on boundary

%% Make rectangular uniform mesh
[x,y] = meshgrid(linspace(a,c,2^(M)),linspace(b,d,2^(M)));
x = reshape(x,[2^(2*M) 1]);
y = reshape(y,[2^(2*M) 1]);
mesh.p = [x,y];
EToV = delaunay(x,y);
mesh.H = EToV;

N=size(EToV,1); % Number of elements
%Nv=size(p,1); % Number of vertex nodes in mesh
Nfaces=size(EToV,2); % Number of faces/element
%VX = p(:,1); % x-coordinates of vertex nodes
%VY = p(:,2); % y-coordinates of vertex nodes

% identify unique vertex nodes
%[p,I,J] = unique(p,'rows');
% change numbering of vertex nodes accordingly to sort order
%EToV = reshape(J(EToV),size(EToV));
% remove duplicate entries in p and update EToV accordingly
%[idx,I,J] = unique(EToV);
%p = p(idx,:);
%EToV = reshape(J,size(EToV));
% Reorder local vertex indexing
%EToV = Reorder(EToV,VX,VY);
h = convhull(mesh.p);
mesh.Nbound = h(1:(end-1));
mesh.BH = [h(1:(end-1)) h(2:end)];

% EToE and EToF
[EToE,EToF]= tiConnect2D(EToV,Nfaces);

% Compute elements on boundary
I = 1:N;
v = (EToE(:,1)' == I | EToE(:,2)' == I | EToE(:,3)' == I);
mesh.Ebound = I(v)';
mesh.Domain = ones(length(EToV),1);
mesh.BDomain = ones(length(mesh.Nbound),1);

% Finally a small plot
figure;
triplot(EToV, x, y)
hold on
for ii = 1:length(mesh.p)
    t = text(mesh.p(ii,1),mesh.p(ii,2),num2str(ii));
    hold on
end

function [EToV] = Reorder(EToV,VX,VY)
% Purpose: Reorder elements to ensure counter-clockwise orientation
ax = VX(EToV(:,1)); ay = VY(EToV(:,1));
bx = VX(EToV(:,2)); by = VY(EToV(:,2));
cx = VX(EToV(:,3)); cy = VY(EToV(:,3));
D = (ax-cx).*(by-cy)-(bx-cx).*(ay-cy);
i = find(D<0);
EToV(i,:) = EToV(i,[1 3 2]);
end
end
