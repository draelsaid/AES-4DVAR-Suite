function [p] = Lopnl(x,Mnl)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L Operator (non-linear) for weak-constraint 4DVAR          %
% [p] - Model errors (4D matrix)                             %
%   x - State variables (4D matrix)                          % 
% Mnl - Function handle containing all propogator info       %
%       (model,scheme,timesteps etc)                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = size(x,2) - 1;   % Number of subwindows determined by size of x
p(:,1) = x(:,1);     % First element of vector p is x_0

for i = 1:s
p(:,i+1) = x(:,i+1) - Mnl(x(:,i));
end

end