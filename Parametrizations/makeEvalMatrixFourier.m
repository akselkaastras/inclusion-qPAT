function evalpar = makeEvalMatrixFourier(xq,maxfreq)
%% Real Fourier basis functions on [0,2*pi]

M = 2*maxfreq+1; % number of basisfunctions
res = length(xq); % number of points

%% When n,m ~= 0
Bn = zeros(res,2*maxfreq);
for n = 1:maxfreq
    index = 2*(n-1);
    Bn(:,index+1) = cos(n*xq);
    Bn(:,index+2) = sin(n*xq);
end

%% Collect
B = [1/sqrt(2*pi)*ones(res,1) 1/sqrt(pi)*Bn];
l = reshape(repmat((1:maxfreq).^2,2,1),1,2*length(1:maxfreq));
Lambda = [0 l]';

%% Save in evalpar struct
evalpar.B = B;
evalpar.Lambda = Lambda;
evalpar.M = M;
evalpar.res = res;
evalpar.maxfreq = maxfreq;
evalpar.xq = xq;