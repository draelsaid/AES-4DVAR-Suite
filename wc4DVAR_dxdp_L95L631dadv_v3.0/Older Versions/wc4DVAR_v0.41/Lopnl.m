function [p] = Lopnl(x,s,n,fm,sch,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L Operator (non-linear) for weak-constraint 4DVAR          %
% [p] - Model errors (4D matrix)                             %
%   x - State variables (4D matrix)                          % 
%   n - Assimilation window time                             %
%  fm - Model function handle                                %
% sch - Numerical scheme function handle                     %
%  dt - timestep                                             %
%   s - number of model errors (s+1 = number of subwindows)  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  p(:,1) = x(:,1); % First element of vector p is x_0
  Mx = x(:,1);

for i=2:(n/s):n+1

   for k=1:(n/s)   % Subwindows
  Mx = Mnl(Mx,fm,sch,dt,n);
   end

  p(:,i) = x(:,i) - Mx;

end

end