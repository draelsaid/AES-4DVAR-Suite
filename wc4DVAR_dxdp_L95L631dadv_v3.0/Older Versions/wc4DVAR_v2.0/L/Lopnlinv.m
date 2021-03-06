function x = Lopnlinv(p,Mnl)
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
s = size(p,2) - 1;   % Number of subwindows determined by size of p
x(:,1) = p(:,1);     % First element of vector p is x_0

for i = 1:s
x(:,i+1) = p(:,i+1) + Mnl(x(:,i));
end

end