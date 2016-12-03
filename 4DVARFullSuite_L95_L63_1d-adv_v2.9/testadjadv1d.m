% Adjoint model test with numerical integration scheme
clear all
close all

addpath( pwd,'L')

% Model settings
N=100; n=120; % Size of spatial and temporal domains
dt = 0.01;    % 0.025 dt => 3 hour time step (L95), 1 time unit is 5 days
dx = 1/N;     % Spatial interval size
C = (-0.8*dt)/(dx)

fnl = @adv1d;
schnl = @lpde;

ftl = @adv1d;
schtl = @lpdetl;

fadj = @adv1d_adj;
schadj = @lpde_adj;

% Initial conditions 
x = randn(N,1); dx = randn(N,1); dy = randn(N,1);

%%%%%%%%%%%%%%%%%%%%%%%% Tests RHS of model equations %%%%%%%%%%%%%%%%%%%%%%
M_dy = ftl(dy,dt); 
Mt_dx = fadj(dx,dt); 

v1 = dot(dx,M_dy); 
v2 = dot(Mt_dx,dy);

adjtest = norm(v1-v2)
disp(sprintf('%5.15f',v1))
disp(sprintf('%5.15f',v2))
% Accuracy should be ~ 10^(-15)

%%%%%%%%%%%%%%%%%%%%%%%% Tests model with propogator (scheme) %%%%%%%%%%%%%%
Mdx = Mnl(dx,fnl,schnl,dt,n);
Mtldx = Mtl(dx,x,fnl,ftl,schtl,schnl,dt,n);
Mtdy = Madj(dy,x,fadj,fnl,schadj,schnl,dt,n);

fwd_prd = dot(dy,Mdx); disp(sprintf('%5.15f',fwd_prd))
adj_prd = dot(Mtdy,dx); disp(sprintf('%5.15f',adj_prd))

fwd_prdtl = dot(dy,Mtldx); disp(sprintf('%5.15f',fwd_prd))
adj_prdtl = dot(Mtdy,dx); disp(sprintf('%5.15f',adj_prd))

Accuracynladj = fwd_prd - adj_prd
Accuracytladj = fwd_prdtl - adj_prdtl
% Accuracy should be ~ 10^(-16)
