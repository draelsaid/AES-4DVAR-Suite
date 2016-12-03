function [x] = lpde(xo,fnl,dt)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % PDE model propagation                               %
   %   x - State variable                                %
   % fnl - Forward model handle                          %
   %  dt - Time-step                                     %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     x  = fnl(xo,dt);
end
