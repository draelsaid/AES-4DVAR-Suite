function [dx] = rk4tl(dx,x,ftl,fnl,dt)  
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Runge-Kutta 4 ODE numerical integration scheme (Tangent-linear) %
   %  dx - Perturbed state                                           %
   %   x - Linearisation state                                       %
   % ftl - Forward tangent linear model handle                       %
   % fnl - Forward Non-linear model handle                           %
   %  dt - Time-step                                                 % 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [~,xmat] = rk4nl(x,fnl,dt); % Call non-linear model to obtain linearisation states 
                              % of each stage of the RK4 scheme. ie. k1,k2,k3,k4 
     zz = dx;
     qq = ftl(xmat(:,1),zz,dt);
     xx = qq;
    
     zz = dx + qq/2;
     qq = ftl(xmat(:,2),zz,dt);
     xx = xx + 2*qq;
     
     zz = dx + qq/2;
     qq = ftl(xmat(:,3),zz,dt);
     xx = xx + 2*qq;
     
     zz = dx + qq;
     qq = ftl(xmat(:,4),zz,dt);
     xx = xx + qq;
     
     dx = dx + xx/6;
end
