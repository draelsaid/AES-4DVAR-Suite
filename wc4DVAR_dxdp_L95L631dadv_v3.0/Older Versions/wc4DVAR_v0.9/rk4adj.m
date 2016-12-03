function dx = rk4adj(dx,x,fadj,fnl,dt)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runge-Kutta 4 ODE numerical integration scheme (Adjoint version)  %
%   dx - Perturbed state variable                                   %
%    x - Linearisation state variable                               %
% fadj - Forward adjoint model handle                               %
%  fnl - Forward Non-linear model handle                            %
%   dt - Time-step                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,xmat] = rk4nl(x,fnl,dt);   % Call non-linear model to obtain linearisation states 
                              % of each stage of the RK4 scheme. ie. k1,k2,k3,k4 
 xx = dx/6;
 
 qq = xx;
 zz = fadj(xmat(:,4),qq,dt);
 dx = dx + zz;
 qq = zz;
  
 qq = qq + 2*xx;
 zz = fadj(xmat(:,3),qq,dt);
 dx = dx + zz;
 qq = 0.5*zz;
     
 qq = qq + 2*xx;
 zz = fadj(xmat(:,2),qq,dt);
 dx = dx + zz;
 qq = 0.5*zz;
     
 qq = qq + xx;
 zz = fadj(xmat(:,1),qq,dt);
 dx = dx + zz;
 
end