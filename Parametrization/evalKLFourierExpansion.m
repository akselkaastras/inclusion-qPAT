function thetaq = evalKLFourierExpansion(cn,evalpar,priorpar) 

% evaluate KL expansion in (xq,yq) using evaluation matrix from evalpar
% and eigenvalues from priorpar
thetaq = evalpar.B * (cn.*(priorpar.lambda).^(1/2));