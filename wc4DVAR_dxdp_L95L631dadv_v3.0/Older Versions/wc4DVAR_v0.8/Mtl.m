function dx = Mtl(dx,x,ftl,fnl,schtl,schnl,dt,n) 
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xtraj(:,1) = x;

 for i=2:n
  xtraj(:,i) = schnl(xtraj(:,i-1),fnl,dt);
 end

 for i=1:n
  dx = schtl(dx,xtraj(:,i),ftl,fnl,dt);
 end

end
