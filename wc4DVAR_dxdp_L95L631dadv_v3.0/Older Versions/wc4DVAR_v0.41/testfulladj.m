% Adjoint model test with numerical integration scheme
clear all
close all

N=40; dt=0.0005; n=200;

J = dt*n

fnl = @l95;
schnl = @rk4nl;

ftl = @l95tl;
schtl = @rk4tl;

fadj = @l95adj;
schadj = @rk4adj;

% Initial conditions 
x = randn(N,1); dx = randn(N,1); dy = randn(N,1);

%%%%%%%%%%%%%%%%%%%%%%%% Tests RHS of model equations %%%%%%%%%%%%%%%%%%%%%%
Mtl_dx = ftl(x,dx,dt); 
Mad_dy = fadj(x,dy,dt); 

v1 = dot(Mtl_dx,dy); 
v2 = dot(Mad_dy,dx);

adjtest = norm(v1-v2)
disp(sprintf('%5.15f',v1))
disp(sprintf('%5.15f',v2))
% Accuracy should be ~ 10^(-16)

%%%%%%%%%%%%%%%%%%%%%%%% Tests model with propogator (scheme) %%%%%%%%%%%%%%
Mdx = Mtl(dx,x,ftl,fnl,schtl,schnl,dt,n);
Mtdy = Madj(dy,x,fadj,fnl,schadj,schnl,dt,n);

fwd_prd = dot(dy,Mdx);
adj_prd = dot(Mtdy,dx);
Accuracy = fwd_prd - adj_prd
% Accuracy should be ~ 10^(-16)