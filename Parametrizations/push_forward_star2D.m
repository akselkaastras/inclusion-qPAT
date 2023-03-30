function gamma = push_forward_star2D(xi_center,theta,xq,yq,priorpar)
% Makes physical parameter consisting of star-shaped inclusions
% out of theta.

% Input: 
%   xi_center [2*ninclusions x 1]: contains centers of inclusions in the
%                                  form [x_1 y_1 x_2 y_2 ... y_ninclusions]
%   xi [priorpar.M x ninclusions]: contains M standard normal coefficients
%                                  for each inclusion
%   priorpar: struct containing information on KL expansion
%   meshpar: struct containing information on FEM mesh

% Reshape xi_center for convenience
xi_center = reshape(xi_center,priorpar.ninclusions,2);

% Initialize
gamma = 0*xq + priorpar.background;

% For each inclusion, we approximate as polygon and add value to FEM-node
% if this is contained in the polygon:
for i = 1:priorpar.ninclusions
    nodes = xi_center(i,:) + [exp(theta(:,i)).*cos(priorpar.angles), exp(theta(:,i)).*sin(priorpar.angles)];
    gamma = gamma + priorpar.v(i)*inpoly2([xq, yq],nodes);
end
