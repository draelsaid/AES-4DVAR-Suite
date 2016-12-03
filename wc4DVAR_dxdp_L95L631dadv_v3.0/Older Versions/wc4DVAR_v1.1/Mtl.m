function [dx,dxtraj] = Mtl(dx,x,ftl,fnl,schtl,schnl,dt,n) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Tangent-linear model propagator                            %
%    dx - Perturbed state variables                           %
%     x - Linearisation state variables                       % 
%  fadj - Linearised Adjoint model function handle            %
%   fnl - Non-linear model function handle                    %
% schtl - Numerical scheme function handle (tl)               %
% schnl - Numerical scheme function handle (nl)               %
%    dt - timestep                                            %
%    n  - numerical integration length                        %
%  traj - denotes 'trajectory' and is a matrix vs time       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xtraj(:,1) = x;

 for i=1:n
  xtraj(:,i+1) = schnl(xtraj(:,i),fnl,dt);
 end

 dxtraj(:,1) = dx;
 
 for i=1:n
  dx = schtl(dx,xtraj(:,i),ftl,fnl,dt);
dxtraj(:,i+1)= dx;
 end

end
