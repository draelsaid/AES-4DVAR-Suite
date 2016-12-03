function [x] = Lopnlinv(p,s,n,fm,sch,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L Operator Inverse (non-linear) for weak-constraint 4DVAR  %
%   [x] - State variables (4D matrix)                        %
%     p - Model errors (4D matrix)                           % 
%     n - Assimilation window time                           %
%    fm - model function handle                              %
%   sch - Numerical scheme function handle                   %
%    dt - timestep                                           %
%     s - number of model errors (s+1 = number of subwindows)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
x(:,1) = p(:,1); % First element of vector p is x_0
Mx = p(:,1);

for i=2:(n/s):n+1
    
  for k=1:(n/s) % Subwindows
 Mx = Mnl(Mx,fm,sch,dt,n);
  end

 x(:,i) = Mx + p(:,i);
 Mx = x(:,i);
end

end