function int = boundint2phi(g,n)
  
  
  % Function calculates the boundary integral of two basis functions from 
  % g(1,:) to g(2,:). 
  %  
  % T. Vilhunen 25.4.2001
  % A. Rasmussen 5.3.2022

% weights and quadrature points
w = [1/2,1/2];
ip = [1/2-1/6*sqrt(3),1/2+1/6*sqrt(3)]; 

% abs quadrature points
x = ip'*(g(2,:)-g(1,:))+g(1,:);
x = x(:,1)+1i*x(:,2);

% compute trigonometric function in quadrature points 
phin = 1/sqrt(2*pi)*exp(1i*double(n)*angle(x));

S = [1-ip(1),1-ip(2);ip(1),ip(2)];

dJt = sqrt((g(2,1)-g(1,1))^2+(g(2,2)-g(1,2))^2); 

int = 1/2*S*phin*dJt;
