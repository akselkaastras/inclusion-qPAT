M = 10000;
x = linspace(0,1,M);
f = @(x) cos(5*x);
df = @(x) abs(5*sin(5*x));
b = cos(5*x+1/2);

m = 1000;
u = 1.*(f(x)>0);
figure(1);
plot(x,u)
eps = zeros(1,M);
eps(33) = 1;
hold on
plot(x,1.*(f(x)+h*eps>0))
hold off
%%
h = 0.1;
dJm = zeros(M,1);
for i = 1:M
    eps = zeros(1,M);
    eps(i) = 1;
    dJm(i) = 1/h * (norm(1.*(f(x)+h*eps>0)-b,2) - norm(1.*(f(x)>0)-b,2));
    i
end
figure(2)
plot(dJm)
hold on
dJmtilde = 0.1*(u-b).*df(x);
plot(dJmtilde)

%%
norm(1.*(f(x)>0)-b,2)
norm(1.*(f(x)-dJmtilde>0)-b,2)
plot(x,1.*(f(x)+dJmtilde>0))