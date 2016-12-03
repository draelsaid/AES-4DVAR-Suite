function [Forward_Product,Adjoint_Product,Accuracy] = test_adjoint(a,L,dx,tassim,dt

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

nx = L/dx;               % Number of grid points
nt = tassim/dt;          % Total time steps

x = randn(nx+1,1);       % Randomly generate state vector
y = randn(nx+1,1);       % Randomly generate state vector
eta = randn(nx+1,nt);    % Randomly generate error state vector

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
end