function dp = LoplT(dx,x,n,fadj,fnl,schadj,schnl,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L Operator Adjoint (linearised) for weak-constraint 4DVAR  %
% [p] - Model errors (4D matrix)                             %
%   x - State variables (4D matrix)                          % 
%   n - Assimilation window time                             %
%  fm - Model function handle                                %
% sch - Numerical scheme function handle                     %
%  dt - timestep                                             %
%   s - number of model errors (s+1 = number of subwindows)  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = size(x,2);     % Number of subwindows determined by size of x
dp(:,1) = dx(:,1); % First element of vector p is x_0

for i = 1:s-1
dp(:,i+1) = dx(:,i+1) - Madj(dx(:,i),x(:,i),fadj,fnl,schadj,schnl,dt,n/(s-1));
end

end