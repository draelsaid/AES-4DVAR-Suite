function dp = Loptl(dx,x,Mtl)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L Operator (linearised) for weak-constraint 4DVAR          %
% [p] - Model errors (4D matrix)                             %
%   x - State variables (4D matrix)                          % 
%   n - Assimilation window time                             %
%  fm - Model function handle                                %
% sch - Numerical scheme function handle                     %
%  dt - timestep                                             %
%   s - number of model errors (s+1 = number of subwindows)  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = size(x,2) - 1;     % Number of subwindows determined by size of x
dp(:,1) = dx(:,1); % First element of vector p is x_0

for i = 1:s
dp(:,i+1) = dx(:,i+1) - Mtl(dx(:,i),x(:,i));
end

end