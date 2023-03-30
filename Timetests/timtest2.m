


r = @(theta) 2*(0.3*sin(theta)+0.1*cos(theta*2));
%r = @(theta) cos(theta);
h = @(theta) exp(r(theta));

atann = @(y,x) 2*atan(y./(sqrt(x.^2+y.^2)+x));

[X,Y] = meshgrid(-1:0.01:1);
xq = X(:);
yq = Y(:);

theta(1) = 0;
theta(2) = 0;
tt = linspace(0,2*pi,1000)';
ht = h(tt);
cost = cos(tt);
sint = sin(tt);

s = [ht.*cost+theta(1) ht.*sint+theta(2)];
figure(1)
plot(s(:,1),s(:,2));
figure(2)
plot(ht);

%% INpolygon is way to slow. We should we really do something better here

for i = 1:100
% arctan or inpolygon algorithm?
% 1. arctan
xc = xq-theta(1);
yc = yq-theta(2);
ind = sqrt(xc.^2+yc.^2) < h(atan2(yc,xc));
val = xc*0;
val(ind) = 1;

figure(1)
%surf(X,Y,reshape(val,[length(X),length(X)]));
%view(2)

% 2. inpolygon
s = [ht.*cost+theta(1) ht.*sint+theta(2)];
ind2 = inpoly2([xq yq],s);

val2 = xc*0;
val2(ind2) = 1;
%figure(2)
%surf(X,Y,reshape(val2,[length(X),length(X)]));
%view(2)
norm(val-val2)
i
end