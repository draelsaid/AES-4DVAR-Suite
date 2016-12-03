%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Purpose: To test adjoint code of adv1d_adj.m
%
%  (c) 2012  University of Reading
%
%  Written  by Adam El-Said
%
%                     y^T A x = x^T A* y
%
% A - model equations
% A* - adjoint model equations (which is just A^T here)
% y,x - random Gaussian generated vectors
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

%% Model parameters
a = 1;                   % Advection model velocity parameter
L = 10;                  % Domain size

dx = 1/20;               % Spatial Resolution
nx = L/dx;               % Number of grid points

%% Assimilation parameters
tassim = 5;              % Assimilation window time length

dt = 1/20;               % Temporal Resolution
nt = tassim/dt;          % Total time steps

C = a*(dt/dx)            % Courant number print

eta=zeros(nx+1,nt);

x = zeros(nx+1,1);
y = zeros(nx+1,1);

% Randomly generated x and y
for i=1:nx+1
    x(i) = randn(1,1);
    y(i) = randn(1,1);
end

% Randomly generated eta at each time step
eta = randn(nx+1,nt);

xo=x;
for n=1:nt-1
[xf]=adv1d(xo,eta(:,n),L,a,1,dt,dx); % 1 time step of model
xo = xf; % x forward (xf)
end

yo=y;
for n=nt-1:-1:1
[xb]=adv1d_adj(yo,eta(:,n),L,a,1,dt,dx); % 1 time step of model
yo = xb; % x backward (xb)
end

Forward_Product = y'*xf
Adjoint_Product = x'*xb

Accuracy = Forward_Product - Adjoint_Product