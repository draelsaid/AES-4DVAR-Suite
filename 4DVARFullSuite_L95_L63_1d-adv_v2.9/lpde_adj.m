function [dx] = lpde_adj(dx,x,fadj,fnl,dt)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Linear PDE model propagation                        %
   %   x  - State variable                               %
   %  dx  - Linearisation state (dummy here)             %
   % fnl  - Forward model handle                         %
   % fadj - Adjoint model handle                         %
   %  dt  - Time-step                                    %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     dx  = fadj(dx,dt);
end
