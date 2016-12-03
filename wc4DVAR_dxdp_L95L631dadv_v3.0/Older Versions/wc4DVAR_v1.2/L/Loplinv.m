function dx = Loplinv(dp,x,n,ftl,fnl,schtl,schnl,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L Operator Inverse (linearised) for weak-constraint 4DVAR  %
%   [x] - State variables (4D matrix)                        %
%     p - Model errors (4D matrix)                           % 
%     n - Assimilation window time                           %
%    fm - model function handle                              %
%   sch - Numerical scheme function handle                   %
%    dt - timestep                                           %
%     s - number of model errors (s+1 = number of subwindows)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = size(dp,2);    % Number of subwindows determined by size of dp
dx(:,1) = dp(:,1); % First element of vector p is x_0

for i = 1:s-1
dx(:,i+1) = dp(:,i+1) + Mtl(dx(:,i),x(:,i),ftl,fnl,schtl,schnl,dt,n/(s-1));
end

end