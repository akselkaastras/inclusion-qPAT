function evalpar = makeEvalMatrixFourier(xq,yq,maxfreq)
%% Real Fourier basis functions on [-1,1]^2



M = (2*maxfreq+1)^2; % number of basisfunctions
res = length(xq); % number of points
xq = pi*xq;
yq = pi*yq;

%% When n,m ~= 0
Bnm = zeros(res,4*maxfreq*maxfreq);
Lambdanm = zeros(4*maxfreq*maxfreq,1);
for n = 1:maxfreq
    for m = 1:maxfreq
        index = (n-1)*4*maxfreq+4*(m-1);
        Bnm(:,index+1) = cos(n*xq).*cos(m*yq);
        Bnm(:,index+2) = cos(n*xq).*sin(m*yq);
        Bnm(:,index+3) = sin(n*xq).*cos(m*yq);
        Bnm(:,index+4) = sin(n*xq).*sin(m*yq);
        Lambdanm(index+(1:4)) = n^2+m^2;
    end
end

%% When n = 0
Bm = zeros(res,2*maxfreq);
Lambdam = zeros(2*maxfreq,1);

for m = 1:maxfreq
    index = 2*(m-1);
    Bm(:,index+1) = cos(m*yq);
    Bm(:,index+2) = sin(m*yq);
    Lambdam(index+(1:2),1) = (m)^2;
end

%% When m = 0
Bn = zeros(res,2*maxfreq);
Lambdan = zeros(maxfreq,1);

for n = 1:maxfreq
    index = 2*(n-1);
    Bn(:,index+1) = cos(n*xq);
    Bn(:,index+2) = sin(n*xq);
    Lambdan(index+(1:2),1) = n^2;
end

%% Collect
B = [1/2*ones(res,1) 1/sqrt(2)*Bm  1/sqrt(2)*Bn Bnm];
Lambda = [0; Lambdan; Lambdam; Lambdanm];

%% Sort in increasing eigenvalues
[Lambda, I] = sort(Lambda);
B = B(:,I);
%% Save in evalpar struct
evalpar.B = B;
evalpar.Lambda = Lambda;
evalpar.M = M;
evalpar.res = res;
evalpar.maxfreq = maxfreq;
evalpar.xq = xq;
evalpar.yq = yq;