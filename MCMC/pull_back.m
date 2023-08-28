function cn = pull_back(gamma,prior_par)

cn = zeros(prior_par.N,2);
% read prior parameters
M = prior_par.M;
N = prior_par.N;
alpha = prior_par.alpha;
tau = prior_par.tau;
q = prior_par.q;

% order coefficients the correct way
odd = 0:N/2-1;
even = -N/2:1:-1;

ei = [odd even];

xi = q*(tau^2+(2*pi)^2*(ei).^2).^(-alpha)';

CN = xi.^(-1/2).*fft(gamma);
cn(:,1) = real(CN);
cn(:,2) = imag(CN);

% reshape to real part and imaginary part in each column
cn = reshape(cn,[2*prior_par.N 1]);