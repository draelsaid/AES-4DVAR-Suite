function [x] = lpdetl(dx,x,ftl,fnl,dt)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % PDE model propagation                               %
   %   x - State variable                                %
   % fnl - Forward model handle                          %
   %  dt - Time-step                                     %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     x  = fnl(dx,dt);
end
