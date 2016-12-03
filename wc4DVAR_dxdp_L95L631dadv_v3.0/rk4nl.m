function [x,xmat] = rk4nl(x,fnl,dt)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Runge-Kutta 4 ODE numerical integration scheme      %
   %   x - State variable                                %
   % fnl - Forward model handle                          %
   %  dt - Time-step                                     %
   %                                                     %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
     zz = x;
     xmat(:,1) = zz;
     qq = fnl(zz,dt);
     xx = qq;
    
     zz = x + qq/2;
     xmat(:,2) = zz;
     qq = fnl(zz,dt);      
     xx = xx + 2*qq;
     
     zz = x + qq/2;
     xmat(:,3) = zz;
     qq = fnl(zz,dt);      
     xx = xx + 2*qq;
     
     zz = x + qq;
     xmat(:,4) = zz;
     qq = fnl(zz,dt);      
     xx = xx + qq;
     
     x = x + xx/6;
end
