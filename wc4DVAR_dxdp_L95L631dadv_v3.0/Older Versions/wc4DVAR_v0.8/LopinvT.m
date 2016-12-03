function [x] = LopinvT(p,s,n,fm,sch,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L Operator Inverse (Transpose) for weak-constraint 4DVAR   %
%   [x] - State variables (4D matrix)                        %
%     p - Model errors (4D matrix)                           % 
%     n - Assimilation window time                           %
%    fm - model function handle                              %
%   sch - Numerical scheme function handle                   %
%    dt - timestep                                           %
%     s - number of model errors (s+1 = number of subwindows)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x(:,n+1) = p(:,n+1); % First element of vector p is x_0
Mx = p(:,n+1);

for i=n:-(n/s):1

  for k=(n/s):-1:1 % Subwindows
 Mx = Mnl(Mx,fm,sch,dt,n);
  end

 x(:,i) = Mx + p(:,i);
 Mx = x(:,i);
end

end