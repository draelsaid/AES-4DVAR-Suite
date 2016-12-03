function dw = LoplinvT(dv,x,n,fadj,fnl,schadj,schnl,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L Operator Adjoint Inverse (linearised) for w-c 4DVAR      %
%   [x] - State variables (4D matrix)                        %
%     p - Model errors (4D matrix)                           % 
%     n - Assimilation window time                           %
%    fm - model function handle                              %
%   sch - Numerical scheme function handle                   %
%    dt - timestep                                           %
%     s - number of model errors (s+1 = number of subwindows)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = size(dv,2) - 1;    % Number of subwindows determined by size of dp
dw(:,s+1) = dv(:,s+1); % First element of vector p is x_0

for i = s:-1:1
dw(:,i) = dv(:,i) + Madj(dw(:,i+1),x(:,i),fadj,fnl,schadj,schnl,dt,n/s);
end

end