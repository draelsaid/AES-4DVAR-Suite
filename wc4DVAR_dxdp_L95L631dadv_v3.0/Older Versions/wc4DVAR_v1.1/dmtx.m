function [Dx,pb] = dmtx(x,sigmab,sigmaq,fcb,fcq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D Matrix Operator for weak-constraint 4DVAR                %
% [Dx] - Model errors (4D matrix)                            %
% [pb] - State variables (4D matrix)                         % 
%   x  - vector                                              %
%  fc  - Correlation function handle                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(x,1); n = size(x,2);

if (nargin == 3) 
  B = sigmab^2*eye(N); Q = sigmaq^2*eye(N);
else
 CB = fcb(N); CQ = fcq(N);    
  B = sigmab^2*CB; Q = sigmaq^2*CQ;
end

Bi = inv(B); Qi = inv(Q);
Brt = sqrt(B); Qrt = sqrt(Q);

% Vector weighting by D
    Dx(:,1) = Bi*x(:,1);
  for i=2:n
    Dx(:,i) = Qi*x(:,i);
  end

% Creation of background and background model errors
    pb(:,1) = Brt*x(:,1);
  for i=2:n
    pb(:,i) = Qrt*x(:,i);
  end
  
end
