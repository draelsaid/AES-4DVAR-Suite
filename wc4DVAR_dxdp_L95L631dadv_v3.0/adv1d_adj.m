function [x]=adv1d_adj(xo,dt)
%% 1D Linear Advection Model (ADJOINT)
%  solved by an upwind scheme with periodic boundary conditions
N = length(xo);        % Spatial domain size determined by user
dx = 1/N;          % Spatial domain interval size
a = -1;              % Wavespeed

uanew = zeros(N,1);    % New solution store
ua=xo;                 % Allocate initial conditions
C = (a*dt)/(dx);       % Courant number

%% Time Loop
for n=1:-1:1

% Periodic Boundary conditions
if a>0
    uanew(N) = (1-C)*(ua(N)) + C*(ua(1));
else
    uanew(1) = (1+C)*(ua(1)) + (-C)*(ua(N));
end

% Upwind Scheme
if a>0
    for i=1:N-1
        uanew(i) = (1-C)*(ua(i)) + C*(ua(i+1));
    end
 else
    for i=2:N
        uanew(i) = (1+C)*(ua(i)) + (-C)*(ua(i-1));
    end
end

% Update
    for i=1:N
     ua(i) = uanew(i);
    end
end
x=ua;
end
