function int = triangint3area(g,mu,order)
  
  % Calculates DA-FEM-matrix integrals phi_i * phi_p * phi_k. 
  % 
  %
  % Modified by Niko H�nninen (June 2017) from version modified by T. Vilhunen 27.9.2001 from test_int2.m (V. Kolehmainen)
  
  % keyboard
  
  if nargin == 2
      order = 3;
  end
  switch order
      case 3
          w = [1/6*ones(3,1)];
          ip = [1/2 0;1/2 1/2;0 1/2];
      case 7
          w = [1/40*ones(3,1); 1/15*ones(3,1);27/120];
          ip = [0 0; 1 0; 0 1; 1/2 0; 1/2 1/2; 0 1/2; 1/3 1/3];
  end
  L = [-1 1 0;-1 0 1];
  Jt = L*g;
  dJt = abs(det(Jt));
  int = 0;
  for ii = 1:length(w)
      S = [1-ip(ii,1)-ip(ii,2);ip(ii,1);ip(ii,2)];
      %   int = int + w(ii)*S*S';
      int1 = w(ii)*mu(1)*S;
      int2 = w(ii)*mu(2)*S;
      int3 = w(ii)*mu(3)*S;
      int = int + int1 + int2 + int3;
      % int = int + w(ii)*S(1)*mu(1)*(S*S');
  end
  int = int*dJt;
