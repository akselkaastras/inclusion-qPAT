function H = smooth_heaviside(x,delta)
H = 0.*(x<=-delta) + 1.*(delta <= x);
H(bitand(-delta<x,x<delta)) = (1/(2*delta).*x(bitand(-delta<x,x<delta))+0.5);
