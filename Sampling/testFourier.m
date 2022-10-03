% test case prior
A = magic(10);
B = ifft2(A);
C = ifft(ifft(A).').';

%% goal;
% Take trigonometric function over [0,1]^2 evaluated in meshgrid([0,1],[0,1]) in NxN
% points. 
% Now make A such that ifft2(A) matches this function.

N = 5;
T = linspace(0,1-1/N,N);
[X,Y] = meshgrid(T,T);

f = @(x,y) 0.*x + exp(2*pi*4i*x).*exp(2*pi*3i*y);
F = f(X,Y);
G2 = 1/N^2*fft2(F);

A = zeros(N,N);
A(2,2) = 1;

G = N^2*ifft2(A);

norm(F-G)
eps = 1e-6;
find(G2<1+eps & G2>1-eps)
G2