function dx = Madj(dx,x,fadj,fnl,schadj,schnl,dt,n) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjoint model propagator                                   %
%    dx - Perturbed state variables                          %
%     x - Linearisation state variables                      % 
%  fadj - Linearised Adjoint model function handle           %
%   fnl - Non-linear model function handle                   %
% schtl - Numerical scheme function handle (adj)             %
% schnl - Numerical scheme function handle (nl)              %
%    dt - timestep                                           %
%    n  - numerical integration length                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xtraj(:,n) = x;

 for i=n:-1:2
  xtraj(:,i-1) = schnl(xtraj(:,i),fnl,dt);
 end

 for i=1:n
  dx = schadj(dx,xtraj(:,i),fadj,fnl,dt);
 end

end
