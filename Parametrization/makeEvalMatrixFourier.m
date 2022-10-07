function evalpar = makeEvalMatrixFourier(xq,yq,M)

M = 2^6; % number of basisfunctions
t = linspace(-1,1,10);
[X,Y] = meshgrid(t);
x = X(:);
y = Y(:);
maxfreq = 4; % maxfreq we loop to
res = length(X(:)); % order 10k

%% When n,m ~= 0
Bnm = zeros(res,4*(maxfreq-1)*(maxfreq-1));
Lambdanm = zeros(4*(maxfreq-1)*(maxfreq-1),1);
for n = 1:maxfreq
    for m = 1:maxfreq
        index = (n-1)*4*maxfreq+4*(m-1);
        Bnm(:,index+1) = cos(n*x).*cos(m*y);
        Bnm(:,index+2) = cos(n*x).*sin(m*y);
        Bnm(:,index+3) = sin(n*x).*cos(m*y);
        Bnm(:,index+4) = sin(n*y).*sin(m*y);
        Lambdanm(index+(1:4)) = n^2+m^2;
    end
end

%% When n = 0
Bm = zeros(res,2*(maxfreq-1));
Lambdam = zeros(maxfreq-1,1);

for m = 2:maxfreq
    index = 2*(m-2);
    Bm(:,index+1) = cos(m*y);
    Bm(:,index+1) = sin(m*y);
    Lambdam(index+(1:2),1) = (m-1)^2;
end

%% When m = 0
Bn = zeros(res,2*(maxfreq-1));
Lambdan = zeros(maxfreq-1,1);

for n = 2:maxfreq
    index = 2*(n-2);
    Bm(:,index+1) = cos(n*x);
    Bm(:,index+1) = sin(n*x);
    Lambdan(index+(1:2),1) = (n-1)^2;
end

%% Collect
B = [ones(res,1) Bm Bn Bnm];
Lambda = [0; Lambdan; Lambdam; Lambdanm];