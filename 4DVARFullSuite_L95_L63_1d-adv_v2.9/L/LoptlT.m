function dw = LoptlT(dv,x,Madj)
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
s = size(x,2) - 1;     % Number of subwindows determined by size of x
dw(:,s+1) = dv(:,s+1);     % First element of vector p is x_0

for i = s:-1:1
dw(:,i) = dv(:,i) - Madj(dv(:,i+1),x(:,i));
end

end