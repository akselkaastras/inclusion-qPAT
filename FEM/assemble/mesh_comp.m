function meshpar = mesh_comp(hmax)

% Create geometry and FEM mesh for the unit disc

% Aksel Rasmussen March 2023

model = createpde();
geometryFromEdges(model,@circleg);
generateMesh(model,"Hmax",hmax);
pdemesh(model); 
axis equal
[p,e,t] = meshToPet(model.Mesh);

% Determine the center points of edges and the corresponding angles
ecenters = (p(:,e(1,:)) + p(:,e(2,:)))/2;
bfii     = angle(ecenters(1,:)+1i*ecenters(2,:));

% Find the set of triangles that have a face on the boundary
% and the index of the corresponding edge.
% We use the fact that each triangle in the mesh has zero, one or two
% vertices on the boundary.
btri_ind = zeros(size(bfii));
for ttt = 1:size(t,2) % Loop over all triangles in the mesh
    xtmp = [];
    % Is vertex 1 of current triangle on the bottom?
    if (p(1,t(1,ttt))^2+p(2,t(1,ttt))^2 > 1-1e-10)
        xtmp = [xtmp, p(:,t(1,ttt))];
    end
    % Is vertex 2 of current triangle on the bottom?
    if (p(1,t(2,ttt))^2+p(2,t(2,ttt))^2 > 1-1e-10)
        xtmp = [xtmp, p(:,t(2,ttt))];
    end
    % Is vertex 3 of current triangle on the bottom?
    if (p(1,t(3,ttt))^2+p(2,t(3,ttt))^2 > 1-1e-10)
        xtmp = [xtmp, p(:,t(3,ttt))];
    end
    % At this point, matrix xtmp has at most 2 columns.
    % If the number of columns is two, we have a triangle that has
    % an edge as one side. Let's determine the edge in question.
    if size(xtmp,2)>1
        tmp_center      = (xtmp(1,1)+xtmp(1,2))/2 + 1i*(xtmp(2,1)+xtmp(2,2))/2;
        tmp_angle       = angle(tmp_center);
        e_ind           = find(abs(bfii-tmp_angle)<1e-8);
        btri_ind(e_ind) = ttt;
    end
end

% Compute lengths of edges
elengths = zeros(size(bfii));
for iii = 1:length(bfii)
    elengths(iii) = norm(p(:,e(1,iii))-p(:,e(2,iii)));
end


% Check the result visually
figure(1)
clf
pdemesh(p,e,t)
axis equal
axis off
%print -dpng mesh.png

% Save result to file
meshpar.p = p;
meshpar.e = e;
meshpar.t = t;

% Centers of mass of triangles
xvec = p(1,:);
trix = mean(xvec(t(1:3,:))); 
yvec = p(2,:);
triy = mean(yvec(t(1:3,:))); 

% angle between points on boundary
z = p(:,e(1,:));
% Angle of nodes
theta = angle(z(1,:) + 1i*z(2,:));
Dfii = diff(sort(theta));
Dfii = Dfii(1);

% Make interpolation points of elements
% Based on EIDORS elem_select
n_interp = 4;
el_dim = 2;
% Get element nodes, and reshape
% need to be of size n_dims_1 x (n_elems*n_dims) for reshape
el_nodes= p(:,t(1:3,:));
el_nodes= reshape(el_nodes, el_dim+1, []);

% Get interpolation matrix
interp= zeros(0,el_dim+1);
for i=0:n_interp
    for j=0:n_interp-i
        interp= [interp;i,j,n_interp-i-j];
    end
end
interp= (interp + 1/(el_dim+1) )/(n_interp+1);

l_interp = size(interp,1);
mdl_pts = interp*el_nodes;
mdl_pts = reshape(mdl_pts, l_interp, length(t(1:3,:)), el_dim);
mdl_pts = permute(mdl_pts, [2,3,1]);

interp_x = squeeze(mdl_pts(:,1,:));
interp_y = squeeze(mdl_pts(:,2,:));


meshpar.trix = trix;
meshpar.triy = triy;
meshpar.btri_ind = btri_ind;
meshpar.Dfii = Dfii;
meshpar.theta = theta;
meshpar.interp_x = interp_x;
meshpar.interp_y = interp_y;
meshpar.hmax = hmax;


