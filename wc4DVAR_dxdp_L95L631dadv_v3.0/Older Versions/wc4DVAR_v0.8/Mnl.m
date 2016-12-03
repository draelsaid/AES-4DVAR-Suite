function xtraj = Mnl(x,fnl,sch,dt,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non-linear model propagator                                %
%    x - State variables                                     % 
%  fnl - Non-linear model function handle                    %
%  sch - Numerical scheme function handle                    %
%   dt - timestep                                            %
%   n  - numerical integration length                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xtraj(:,1) = x;

for i=2:n+1
xtraj(:,i) = sch(xtraj(:,i-1),fnl,dt);
end


end
