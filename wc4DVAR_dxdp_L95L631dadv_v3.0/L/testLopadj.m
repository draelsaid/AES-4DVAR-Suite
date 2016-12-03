% Adjoint model test with numerical integration scheme
clear all
close all

addpath(pwd,'L')

N=50; dt=0.02; n=200;

fnl = @adv1d;       schnl = @lpde;
ftl = @adv1d;       schtl = @lpdetl;
fadj= @adv1d_adj;  schadj = @lpde_adj;

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

result1 = abs(fwd_prd - adj_prd); disp(sprintf('%5.15f',result1))
result1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Linvdp = Loptlinv(dx,x,ftlm);
Linvtdyp = LoptlinvT(dy,x,fadjm);

fwd_prd2 = sum(sum(dy.*Linvdp)); disp(sprintf('%5.15f',fwd_prd2))
adj_prd2 = sum(sum(Linvtdyp.*dx)); disp(sprintf('%5.15f',adj_prd2))

result2 = abs(fwd_prd2 - adj_prd2); disp(sprintf('%5.15f',result2))
result2

