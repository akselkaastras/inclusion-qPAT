function [xinterp,w,gamma] = push_forward1D(cn,prior_par)

% read prior parameters
v = prior_par.v; % values at each interface (M)
c = prior_par.c; % contour levels (M + 1)
M = prior_par.M; % number of interfaces
N = prior_par.N; % number of fourier coefficients
Ninterp = prior_par.Ninterp;
xinterp = prior_par.xinterp;

% Matern covariance parameters
alpha = prior_par.alpha;
tau = prior_par.tau;
q = prior_par.q;

% reshape to real part and imaginary part in each column
cn = reshape(cn,[prior_par.N 1]);

% order coefficients the correct way
odd = 0:N/2-1;
even = -N/2:1:-1; 

ei = [odd even];

% decaying weights
xi = q*(tau^2+(2*pi)^2*(ei).^2).^(-alpha)';

% Fourier coefficients
a = sqrt(xi).*cn;

% Append zeros for interpolation
nyqst = ceil((N+1)/2);
b = [a(1:nyqst) ; zeros(Ninterp-N,1) ; a(nyqst+1:N)];
b(nyqst) = b(nyqst)/2;
b(nyqst+Ninterp-N) = b(nyqst);

% Function in periodic basis (1,cos(x),sin(x),cos(2x),sin(2x),...
w = ifft(b);
w = real(w) + imag(w);

w = w * Ninterp / N;
w = [w; w(1)]; % append first value as last


% identify each levelset and insert inclusion
gamma = 0*xinterp;
for i = 1:M
    gamma = gamma + v(i)*bitand(c(i)<w,w<=c(i+1))';
end