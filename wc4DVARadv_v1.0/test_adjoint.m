function [Forward_Product,Adjoint_Product,Accuracy] = test_adjoint(a,nx,dx,nt,dt)

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
% y,x - random state vectors
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = randn(nx,1);       % Generate random state vector
y = randn(nx,1);       % Generate random state vector

xo=x;
for n=1:nt
[xf]=adv1d(xo,a,nx,dx,1,dt); % 1 time step of model
xo = xf; % x forward (xf)
end

yo=y;
for n=nt:-1:1
[xb]=adv1d_adj(yo,a,nx,dx,1,dt); % 1 time step of model
yo = xb; % x backward (xb)
end

Forward_Product = y'*xf; Adjoint_Product = x'*xb; Accuracy = Forward_Product - Adjoint_Product;

end