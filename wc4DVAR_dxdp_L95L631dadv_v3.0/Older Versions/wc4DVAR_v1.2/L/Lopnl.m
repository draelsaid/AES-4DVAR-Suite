function p = Lopnl(x,n,fm,sch,dt)
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
s = size(x,2);   % Number of subwindows determined by size of x
p(:,1) = x(:,1); % First element of vector p is x_0

for i = 1:s-1
p(:,i+1) = x(:,i+1) - Mnl(x(:,i),fm,sch,dt,n/(s-1));
end

end