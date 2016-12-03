function [x]=adv1d(xo,a,nx,dx,tstep,dt)
%% 1D Linear Advection Model (IMPERFECT MODEL)
%  solved by an upwind scheme with periodic boundary condition

unew = zeros(nx,1);     % New solution store
u=xo;                   % Allocate initial conditions
C = (a*dt)/(dx);        % Courant number

%% Time Loop
for n=1:tstep

% Periodic Boundary conditions
if a>0
    unew(1) = (1-C)*(u(1)) + (C)*(u(nx));
else
    unew(nx) = (1+C)*(u(nx)) + (-C)*(u(1));
end

% Upwind Scheme
if a>0
    for i=2:nx
        unew(i) = C*(u(i-1)) + (1-C)*(u(i));
    end
 else
     for i=1:nx-1
        unew(i) = (1+C)*(u(i)) + (-C)*(u(i+1));
     end
end
% Update
    for i=1:nx
     u(i) = unew(i);
    end
end
x=u;
end
