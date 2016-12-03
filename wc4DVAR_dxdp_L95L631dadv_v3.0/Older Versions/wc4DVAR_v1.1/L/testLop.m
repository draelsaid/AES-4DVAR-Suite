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

%%%%%%%%%%%%%%%%%%%%%%%% (INVERSE TEST) Tests non-linear L %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial conditions 
x = randn(N,s+1);

p = Lopnl(x,n,fnl,schnl,dt);
x1 = Lopnlinv(p,n,fnl,schnl,dt);

t = x-x1; % Should =0

%%%%%%%%%%%%%%%%%%%%%%%% (INVERSE TEST) Tests linearised L %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial conditions 
x = randn(N,s+1); dx = randn(N,s+1);

dp = Lopl(dx,x,n,ftl,fnl,schtl,schnl,dt);
dx1 = Loplinv(dp,x,n,ftl,fnl,schtl,schnl,dt);

t2 = dx-dx1; % Should =0

%%%%%%%%%%%%%%%%%%%%%%%% (INVERSE TEST) Tests Adjoint linearised L %%%%%%%%%%%%%%%%%%%%%%

% Initial conditions 
x = randn(N,s+1); dx = randn(N,s+1);

dp = LoplT(dx,x,n,fadj,fnl,schadj,schnl,dt);
dx2 = LoplinvT(dp,x,n,fadj,fnl,schadj,schnl,dt);

t3 = dx-dx2; % Should =0

