function xtraj = Mnlsc(x,Mnl,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non-linear model propagator for SC4DVAR                    %
%    x - State variables                                     % 
%  fnl - Non-linear model function handle                    %
%  sch - Numerical scheme function handle                    %
%   dt - timestep                                            %
%   n  - numerical integration length                        %
%  traj - denotes 'trajectory' and is a matrix vs time       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xtraj(:,1) = x;

for i=1:n    
xtraj(:,i+1) = Mnl(xtraj(:,i));
end

end