% L operator test
clear all
close all

addpath( pwd,'L')

N=40; dt=0.025; n=200; s=n/2;

fnl = @l95;
schnl = @rk4nl;

ftl = @l95tl;
schtl = @rk4tl;

fadj = @l95adj;
schadj = @rk4adj;

fnlm = @(x) Mnl(x,fnl,schnl,dt,1);
ftlm = @(x,dx) Mtl(dx,x,ftl,fnl,schtl,schnl,dt,1); 
fadjm = @(y,dy) Madj(dy,y,fadj,fnl,schadj,schnl,dt,1);

%%%%%%%%%%%%%%%%%%%%%%%% (INVERSE TEST) Tests non-linear L %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial conditions 
x = randn(N,s+1);

p = Lopnl(x,fnlm);
x1 = Lopnlinv(p,fnlm);

t = x-x1; % Should =0

%%%%%%%%%%%%%%%%%%%%%%%% (INVERSE TEST) Tests linearised L %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial conditions 
x = randn(N,s+1); dx = randn(N,s+1);

dp = Loptl(dx,x,ftlm);
dx1 = Loptlinv(dp,x,ftlm);

t2 = dx-dx1; % Should =0

%%%%%%%%%%%%%%%%%%%%%%%% (INVERSE TEST) Tests Adjoint linearised L %%%%%%%%%%%%%%%%%%%%%%

% Initial conditions 
x = randn(N,s+1); dx = randn(N,s+1);

dp = LoptlT(dx,x,fadjm);
dx2 = LoptlinvT(dp,x,fadjm);

t3 = dx-dx2; % Should =0

