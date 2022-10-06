function gamma = push_forward_levelset2D(theta,priorpar)

% read prior parameters
%v = priorpar.v; % values at each interface (M)
%c = priorpar.c; % contour levels (M + 1)
%M = priorpar.M; % number of interfaces


% identify each levelset and insert inclusion
gamma = 0*theta + 1;
for i = 1:priorpar.M
    gamma = gamma + priorpar.v(i)*bitand(priorpar.c(i)<theta,theta<=priorpar.c(i+1));
end

