% Kite and cushion parametrizations
t = linspace(0,2*pi,1000);
kitecurve = [cos(t)+0.65*cos(2*t)-0.65; 1.5*sin(t)];
r = sqrt(0.8+0.8*(cos(4*t)-1).^2);
cushioncurve = [r.*cos(t); r.*sin(t)]; 
curves = [0.18*kitecurve+[0.4;-0.4]; 0.12*cushioncurve+[-0.4;0.4]];
values = [0.2,0.4,0.1];
noiselevel = 0.01;
seed = 0;

datapar = make_data(curves,values,noiselevel,fine_hmax,meshpar,seed); 