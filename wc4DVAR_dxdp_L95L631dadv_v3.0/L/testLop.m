% L operator test
clear all
close all

addpath( pwd,'L')

N=50; dt=0.02; n=200;

fnl = @adv1d;
schnl = @lpde;

ftl = @adv1d;
schtl = @lpdetl;

fadj = @adv1d_adj;
schadj = @lpde_adj;

fnlm = @(x) Mnl(x,fnl,schnl,dt,1);
ftlm = @(x,dx) Mtl(dx,x,ftl,fnl,schtl,schnl,dt,1); 
fadjm = @(y,dy) Madj(dy,y,fadj,fnl,schadj,schnl,dt,1);

%%%%%%%%%%%%%%%%%%%%%%%% (INVERSE TEST) Tests non-linear L %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial conditions 
x = randn(N,n+1);

p = Lopnl(x,fnlm);
x1 = Lopnlinv(p,fnlm);

t1 = x-x1; % Should =0
normt1 = norm(t1)

p = randn(N,n+1);

x = Lopnlinv(p,fnlm);
p1 = Lopnl(x,fnlm);

t1b = p-p1; % Should =0
normt1b = norm(t1b)

%%%%%%%%%%%%%%%%%%%%%%%% (INVERSE TEST) Tests linearised L %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial conditions 
x = randn(N,n+1); dx = randn(N,n+1);

dp = Loptl(dx,x,ftlm);
dx1 = Loptlinv(dp,x,ftlm);

t2 = dx-dx1;
normt2 = norm(t2) % Should =0

p = randn(N,n+1); dp = randn(N,n+1);

dx = Loptlinv(dp,p,ftlm);
dp1 = Loptl(dx,p,ftlm);

t2b = dp - dp1; 
normt2b = norm(t2b) % Should =0

%%%%%%%%%%%%%%%%%%%%%%%% (INVERSE TEST) Tests Adjoint linearised L %%%%%%%%%%%%%%%%%%%%%%

% Initial conditions 
x = randn(N,n+1); dx = randn(N,n+1);

dp = LoptlT(dx,x,fadjm);
dx2 = LoptlinvT(dp,x,fadjm);

t3 = dx-dx2; % Should =0
normt3 = norm(t3)

p = randn(N,n+1); dp = randn(N,n+1);

dx = LoptlinvT(dp,p,fadjm);
dp1 = LoptlT(dx,p,fadjm);

t3b = dp - dp1; % Should =0
normt3b = norm(t3b)
