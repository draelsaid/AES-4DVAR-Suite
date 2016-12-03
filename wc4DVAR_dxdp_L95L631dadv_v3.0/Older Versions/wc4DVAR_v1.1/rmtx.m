function [Rx,y] = rmtx(x,sigmar,fcr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R Matrix Operator                                          %
% [Rx] - Model errors (4D matrix)                            %
% [pb] - State variables (4D matrix)                         % 
%   x  - vector                                              %
%  fc  - Correlation function handle                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(x,1); n = size(x,2);

if (nargin == 2) 
  R = sigmar^2*eye(N);
else
 CR = fcr(N);    
  R = sigmar^2*CR;
end

Ri = inv(R); Rrt = sqrt(R);

% Vector weighting by R
  for i=1:n
    Rx(:,i) = Ri*x(:,i);
  end

% Creation of observation noise
  for i=1:n
    y(:,i) = Rrt*x(:,i);
  end
  
end
