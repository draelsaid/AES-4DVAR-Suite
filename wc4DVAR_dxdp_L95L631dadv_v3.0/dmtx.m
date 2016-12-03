function [Dix,pb,B,Q] = dmtx(x,sigmab,sigmaq,fcb,fcq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D Matrix Operator for weak-constraint 4DVAR                %
% [Dx] - Model errors (4D matrix)                            %
% [pb] - State variables, with model errors (4D matrix)      % 
%   x  - vector                                              %
%  fcb  - Correlation function, function handle (background) %
%  fcq  - Correlation function, function handle (model error)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(x,1); n = size(x,2);

if (nargin == 3)
  B = sigmab^2*eye(N); Q = sigmaq^2*eye(N);
else
    
atypeb = iterchk(fcb);    
atypeq = iterchk(fcq);

if (isempty(fcb))
  B = sigmab^2*eye(N);
end

if (isempty(fcq))
  Q = sigmaq^2*eye(N);
end

if strcmp(atypeb,'matrix')    
 CB = fcb; 
  B = sigmab^2*CB;
end

if strcmp(atypeq,'matrix')
 CQ = fcq;    
  Q = sigmaq^2*CQ;
end

if strcmp(atypeb,'function')
 CB = fcb(N); 
  B = sigmab^2*CB;  
end

if strcmp(atypeq,'function')
 CQ = fcq(N); 
  Q = sigmaq^2*CQ;
end

end
Bi = inv(B); Qi = inv(Q);
Brt = sqrtm(B); Qrt = sqrtm(Q);

% Vector weighting by D
    Dix(:,1) = Bi*x(:,1);
  for i=2:n
    Dix(:,i) = Qi*x(:,i);
  end

% Creation of background and background model errors
    pb(:,1) = Brt*x(:,1);
  for i=2:n
    pb(:,i) = Qrt*x(:,i);
  end
  
end
