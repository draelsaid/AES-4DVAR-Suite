% Testing the tangent-linear model with the numerical integration scheme
clear all
close all

addpath( pwd,'L')

N=40; dt=0.025; n=200; 

fnl = @l95;
schnl = @rk4nl;

ftl = @l95tl;
schtl = @rk4tl;

fLnl = @(x)    Lopnl(x,n,fnl,schnl,dt);
fLtl = @(dx,x) Lopl(dx,x,n,ftl,fnl,schtl,schnl,dt);

% Correctness Test for Lnl and Ltl
for i=1:16
alpha = 10^(1-i);
dx = alpha*randn(N,n+1);
x = randn(N,n+1);
xdx = x + dx;

Ldp = fLtl(dx,x);
lpdp = fLnl(xdx);
lp   = fLnl(x);

alphadx(i)  = norm(lpdp - lp)/norm(Ldp);
alphadx2(i) = alphadx(i)-1;
alphadx3(i) = norm(lpdp - lp - Ldp);
end

semilogy(abs(alphadx)); figure
semilogy(abs(alphadx2)); figure 
semilogy(abs(alphadx3))

% Correctness Test for Lnlinv and Ltlinv
fLnl = @(p)    Lopnlinv(p,n,fnl,schnl,dt);
fLtl = @(dp,p) Loplinv(dp,p,n,ftl,fnl,schtl,schnl,dt);

for i=1:16
alpha = 10^(1-i);
dp = alpha*randn(N,n+1);
p = randn(N,n+1);
pdp = p + dp;


lxdx = fLnl(pdp);
lx   = fLnl(p);
Ldx = fLtl(dp,lx);
 
alphadx(i)  = norm(lxdx - lx)/norm(Ldx);
alphadx2(i) = alphadx(i)-1;
alphadx3(i) = norm(lxdx - lx - Ldx);
end

semilogy(abs(alphadx)); figure
semilogy(abs(alphadx2)); figure 
semilogy(abs(alphadx3))
