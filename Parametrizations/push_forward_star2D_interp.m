function gamma = push_forward_star2D_interp(xi_center,theta,xq,yq,n,HN,priorpar)
% Makes physical parameter consisting of star-shaped inclusions
% out of theta.


% Reshape xi_center for convenience
xi_center = reshape(xi_center,priorpar.ninclusions,2);

% Initialize
gamma = zeros(HN,1) + priorpar.background;

% For each inclusion, we approximate as polygon and add value to FEM-node
% if this is contained in the polygon:
for i = 1:priorpar.ninclusions
    nodes = xi_center(i,:) + [exp(theta(:,i)).*cos(priorpar.angles), exp(theta(:,i)).*sin(priorpar.angles)];
    
    gamma = gamma + priorpar.v(i)*mean(reshape(inpoly2([xq, yq],nodes),[HN n]),2);
end
