function gamma = push_forward_levelset2D_smooth(theta,priorpar)

delta = priorpar.delta;

% identify each levelset and insert inclusion
gamma = 0*theta;
for i = 1:priorpar.ninterface
    gamma = gamma + priorpar.v(i)*(smooth_heaviside(theta-priorpar.c(i),delta)-smooth_heaviside(theta-priorpar.c(i+1),delta));
end

