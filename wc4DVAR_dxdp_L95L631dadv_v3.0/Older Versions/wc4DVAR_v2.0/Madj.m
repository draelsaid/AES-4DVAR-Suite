function [dx,dxtraj] = Madj(dx,x,fadj,fnl,schadj,schnl,dt,n)
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
%  traj - denotes 'trajectory' and is a matrix vs time       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xtraj(:,n+1) = x;

 for i=n+1:-1:2
  xtraj(:,i-1) = schnl(xtraj(:,i),fnl,dt); % Non-linear model run on linearisation point
 end

 dxtraj(:,1) = dx;
 
 for i=1:n
  dx = schadj(dx,xtraj(:,i+1),fadj,fnl,dt);
  dxtraj(:,i+1)=dx;
 end

end
