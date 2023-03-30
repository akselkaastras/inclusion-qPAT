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

f = @(x,y) 0.*x + exp(-1*1i*2*pi*x).*exp(2*1i*2*pi*y);
F = f(X,Y);
G2 = 1/N^2*fft2(F)
fftshift(G2)

% eigenvector is 
%% 

A = zeros(N,N);
A(2,2) = 1;
G = N^2*ifft2(A);

% G corresponds to the signal exp(1i*2*pi*2*x)*exp(1i*2*pi*2*y)

%% Fill in correctly.
k = 2;
c = 4;
M = 2^k+1;
N = 2^c+1;
Lambda = zeros(N,N);
% Eigenvalues of q*(tau^2 - Delta)^{-alpha} is q*(tau^2 + lambda)^{-alpha}
maxfreq = (M-1)/2;
[X,Y] = meshgrid(-maxfreq:maxfreq);
A = X.^2+Y.^2;
hm
%%
Lambda(-maxfreq:maxfreq) = A;
%Lambda
ifftshift(A)

%% Identify where zeros go from interpolation code