[X,Y] = meshgrid(-1:0.01:1);
xq = reshape(X,[numel(X) 1]);
yq = reshape(Y,[numel(Y) 1]);
res = length(xq);
maxfreq = 3;
M = (2*maxfreq+1)^2;
f = fourierBasis(maxfreq,2);
evalpar = makeEvalMatrixFourier(xq,yq,maxfreq);
lambda = (1:maxfreq).^2;
lambda = [0, reshape([lambda; lambda], [1 2*numel(lambda)])];
lambda = reshape(lambda' + lambda,[1 numel(lambda).^2]);
[lambda, I] = sort(lambda);

Bnm = zeros(res,M-1);
for j = 1:length(xq)
    Bnm(j,:) = f(xq(j),yq(j));
end
Bnm = [0.5*ones(res,1) Bnm];
Bnm = Bnm(:,I);