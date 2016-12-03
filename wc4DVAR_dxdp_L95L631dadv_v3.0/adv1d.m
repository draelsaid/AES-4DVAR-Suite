function [x]=adv1d(xo,dt)
%% 1D Linear Advection Model
%  solved using upwind scheme with periodic boundary conditions
N = length(xo);     % Spatial domain size determined by user
dx = 1/N;           % Spatial domain interval size
a = -1;             % Wavespeed

unew = zeros(N,1);      % New solution store
u = xo;                 % Allocate initial conditions
C = (a*dt)/(dx);        % Courant number

%% Time Loop
for n=1:1

% Periodic Boundary conditions
if a>0
    unew(1) = (1-C)*(u(1)) + (C)*(u(N));
else
    unew(N) = (1+C)*(u(N)) + (-C)*(u(1));
end

% Upwind Scheme
if a>0
    for i=2:N
        unew(i) = C*(u(i-1)) + (1-C)*(u(i));
    end
 else
     for i=1:N-1
        unew(i) = (1+C)*(u(i)) + (-C)*(u(i+1));
     end
end

% Update
    for i=1:N
     u(i) = unew(i);
    end
end
x = u;
end
