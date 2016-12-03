% Adjoint model test with numerical integration scheme
clear all
close all

N=40; dt=0.025; n=200;

fnl = @l95;
schnl = @rk4nl;

ftl = @l95tl;
schtl = @rk4tl;

fadj = @l95adj;
schadj = @rk4adj;

% Initial conditions 
x = randn(N*(n+1),1); dx = randn(N*(n+1),1); dp = randn(N*(n+1),1);

%%%%%%%%%%%%%%%%%%%%%%%% Tests model with propogator (scheme) %%%%%%%%%%%%%%
p = Lop(x,s,n,fnl,rk4nl,dt);
[p] = LopT(x,s,n,fadj,rk4adj,dt);

fwd_prd = dot(dy,Mdx);
adj_prd = dot(Mtdy,dx);
Accuracy = fwd_prd - adj_prd
% Accuracy should be ~ 10^(-16)