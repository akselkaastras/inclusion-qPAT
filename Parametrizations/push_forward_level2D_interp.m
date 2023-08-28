function gamma = push_forward_level2D_interp(theta,n,HN,priorpar)
% Makes physical parameter


% Initialize
gamma = zeros(HN,1);

% For each inclusion, we approximate as polygon and add value to FEM-node
% if this is contained in the polygon:
for i = 1:priorpar.ninterface
    
    gamma = gamma + priorpar.v(i)*mean(reshape(bitand(priorpar.c(i)<theta,theta<=priorpar.c(i+1)),[HN n]),2);
end

