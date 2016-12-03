% Adjoint model test with numerical integration scheme
clear all
close all

addpath(pwd,'L')

N=40; dt=0.025; n=200;

fnl = @l95;      schnl = @rk4nl;
ftl = @l95tl;    schtl = @rk4tl;
fadj= @l95adj;  schadj = @rk4adj;

fnlm = @(x) Mnl(x,fnl,schnl,dt,1);
ftlm = @(dx,x) Mtl(dx,x,ftl,fnl,schtl,schnl,dt,1);
fadjm = @(dy,y) Madj(dy,y,fadj,fnl,schadj,schnl,dt,1);

% Initial conditions 
x = randn(N,n+1); dx = randn(N,n+1); dy = randn(N,n+1);

%%%%%%%%%%%%%%%%%%%%%%%% Tests model with propogator (scheme) %%%%%%%%%%%%%%
Ldx = Loptl(dx,x,ftlm);
Ltdy = LoptlT(dy,x,fadjm);

fwd_prd = sum(sum(dy.*Ldx)); disp(sprintf('%5.15f',fwd_prd))
adj_prd = sum(sum(Ltdy.*dx)); disp(sprintf('%5.15f',adj_prd))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ldp = Loptlinv(dx,x,ftlm);
Ltdyp = LoptlinvT(dy,x,fadjm);

fwd_prd = sum(sum(dy.*Ldp)); disp(sprintf('%5.15f',fwd_prd))
adj_prd = sum(sum(Ltdyp.*dx)); disp(sprintf('%5.15f',adj_prd))

