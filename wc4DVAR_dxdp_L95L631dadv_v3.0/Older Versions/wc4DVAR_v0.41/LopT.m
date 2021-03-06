function [dp] = LopT(x,dx,s,n,fnl,fadj,schnl,schadj,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L Operator (Tranpose) for weak-constraint 4DVAR            %
% [p]    - Model errors (4D matrix)                          %
%   x    - State variables (4D matrix)                       %
%   n    - Assimilation window time                          %
% fnl    - Non-linear model function handle                  %
% fadj   - Adjoint model function handle                     %
% schnl  - Non-linear model numerical scheme function handle %
% schadj - Adjoint model numerical scheme function handle    %
%  dt    - Timestep                                          %
%   s    - Number of model errors (s+1 = subwindows)         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p(:,n+1) = x(:,n+1); % First element of vector p is x_0
Mx = x(:,n+1);

for i=n+1:-(n/s):1

   for k=(n/s):1     % Subwindows
  Mx = Madj(dx,x,fadj,fnl,schadj,schnl,dt,n) ;
   end

p(:,i) = x(:,i) - Mx;

end

end