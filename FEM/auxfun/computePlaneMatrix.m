function G = computePlaneMatrix(x,y)

% G is the matrix such that given z = [z1,z2,z3]
% G * z = v = [a,b,c], where z = ax + by + c

G = [x, y, ones(size(x))];
G = inv(G);