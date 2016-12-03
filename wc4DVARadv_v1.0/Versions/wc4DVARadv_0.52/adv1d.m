function [x]=adv1d(xo,eta,nx,a,tstep,dt,dx)
%% 1D Linear Advection Model (IMPERFECT MODEL)
%  solved by an upwind scheme with periodic boundary condition

unew = zeros(nx+1,1);   % New solution store
u=xo;                   % Allocate initial conditions
C = (a*dt)/(dx);        % Courant number

%% Time Loop
for n=1:tstep

% Periodic Boundary conditions
if a>0
    unew(1) = (1-C)*(u(1)) + C*(u(nx+1)) + eta(1);
else
    unew(nx+1) = (1+C)*(u(nx+1)) + (-C)*(u(1)) + eta(nx+1);
end

% Upwind Scheme
if a>0
    for i=2:nx+1
        unew(i) = C*(u(i-1)) + (1-C)*(u(i)) + eta(i);
    end
 else
     for i=1:nx
        unew(i) = (1+C)*(u(i)) + (-C)*(u(i+1)) + eta(i);
     end
 end

% Update
    for i=1:nx+1
     u(i) = unew(i);
    end

end

x=u;

end
