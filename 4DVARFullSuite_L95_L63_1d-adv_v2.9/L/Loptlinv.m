function dx = Loptlinv(dp,x,Mtl)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L Operator Inverse (linearised) for weak-constraint 4DVAR  %
%   [dx] - State variables (4D matrix)                        %
%     p - Model errors (4D matrix)                           % 
%     n - Assimilation window time                           %
%    fm - model function handle                              %
%   sch - Numerical scheme function handle                   %
%    dt - timestep                                           %
%     s - number of model errors (s+1 = number of subwindows)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = size(dp,2) - 1;    % Number of subwindows determined by size of dp
dx(:,1) = dp(:,1); % First element of vector p is x_0

for i = 1:s
dx(:,i+1) = dp(:,i+1) + Mtl(dx(:,i),x(:,i));
end

end