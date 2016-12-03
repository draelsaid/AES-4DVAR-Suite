% L operator test
clear all
close all

addpath( pwd,'L')

N=40; dt=0.025; n=10; s=n/2;

fnl = @l95;
schnl = @rk4nl;

ftl = @l95tl;
schtl = @rk4tl;

fadj = @l95adj;
schadj = @rk4adj;

% Initial conditions 
x = randn(N,s);

%%%%%%%%%%%%%%%%%%%%%%%% Tests non-linear L %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = Lopnl(x,n,fnl,schnl,dt);
x1 = Lopnlinv(p,n,fnl,schnl,dt);

t = x-x1; % Should =0

%%%%%%%%%%%%%%%%%%%%%%%% Tests linearised L %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial conditions 
x = randn(N,s); dx = randn(N,s);

dp = Lopl(dx,x,n,ftl,fnl,schtl,schnl,dt);
dx1 = Loplinv(dp,x,n,ftl,fnl,schtl,schnl,dt);

t2 = dx-dx1; % Should =0

%%%%%%%%%%%%%%%%%%%%%%%% Tests Adjoint linearised L %%%%%%%%%%%%%%%%%%%%%%

% Initial conditions 
x = randn(N,s); dx = randn(N,s);

dp = LoplT(dx,x,n,fadj,fnl,schadj,schnl,dt);
dx2 = LoplinvT(dp,x,n,fadj,fnl,schadj,schnl,dt);

t3 = dx-dx2; % Should =0

%%%%%%%%%%%%%%%%%%%%%%%% Tests if L^T is adjoint of L %%%%%%%%%%%%%%%%%%%%

% Initial conditions 
x = randn(N,s); dxo = randn(N,s);

dp = LoplT(dxo,x,n,fadj,fnl,schadj,schnl,dt);
dx = LoplinvT(dp,x,n,fadj,fnl,schadj,schnl,dt);

dp = Lopl(dx,x,n,ftl,fnl,schtl,schnl,dt);
dx3 = Loplinv(dp,x,n,ftl,fnl,schtl,schnl,dt);

t4 = dxo-dx3; % Should = 0
