clear all
close all

N=40; 

% Initial conditions
dx = randn(N,1); 
dy = randn(N,1);
x = randn(N,1); 

[Mtl_dx] = l95tl(x,dx); 
[Mad_dy] = l95adj(x,dy);

v1 = dot(Mtl_dx,dy); 
v2 = dot(Mad_dy,dx);

adjtest = norm(v1-v2)
disp(sprintf('%5.15f',v1))
disp(sprintf('%5.15f',v2))