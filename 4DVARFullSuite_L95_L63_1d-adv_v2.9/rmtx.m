function [Rx,y] = rmtx(x,sigmao,fcr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R Matrix Operator                                          %
% [Rx] - vector x weighted by R                              %
% [y]  - Obs variables (4D matrix)                           % 
%   x  - vector                                              %
%  fcr - Correlation function, function handle               %
%sigmao- Obs error variance                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(x,1); n = size(x,2);

if (nargin == 2) 
  R = sigmao^2*eye(N);
else
 CR = fcr(N);    
  R = sigmao^2*CR;
end

Ri = inv(R); Rrt = sqrtm(R);

% Vector weighting by R
  for i=1:n
    Rx(:,i) = Ri*x(:,i);
  end

% Creation of observation noise
  for i=1:n
    y(:,i) = Rrt*x(:,i);
  end
  
end
