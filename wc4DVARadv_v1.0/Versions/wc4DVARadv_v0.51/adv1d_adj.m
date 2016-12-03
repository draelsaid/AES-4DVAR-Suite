function [mx]=adv1d_adj(xo,eta,nx,a,tstep,dt,dx)

%% 1D Linear Advection Model (ADJOINT OF IMPERFECT MODEL)
%  solved by an upwind scheme with periodic boundary condition

uanew = zeros(nx+1,1);   % New solution store
ua=xo;                   % Allocate initial conditions
C = (a*dt)/(dx);         % Courant number

%% Time Loop
for n=tstep:-1:1

% Periodic Boundary conditions
if a>0
    uanew(nx+1) = (1-C)*(ua(nx+1)) + C*(ua(1)) + eta(nx+1);
else
    uanew(1) = (1+C)*(ua(1)) + (-C)*(ua(nx+1)) + eta(1);
end

% Upwind Scheme
if a>0
    for i=1:nx
        uanew(i) = (1-C)*(ua(i)) + C*(ua(i+1)) + eta(i);
    end
 else
    for i=2:nx+1
        uanew(i) = (1+C)*(ua(i)) + (-C)*(ua(i-1)) + eta(i);
    end
 end

% Update
    for i=1:nx+1
     ua(i) = uanew(i);
    end

end

mx=ua;

end
