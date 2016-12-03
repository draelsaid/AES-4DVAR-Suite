% Testing the tangent-linear model with the numerical integration scheme
clear all
close all

addpath( pwd,'L')

N=40; dt=0.025; n=20; 

fnl = @l95; schnl = @rk4nl;
ftl = @l95tl; schtl = @rk4tl;

fnlm = @(x) Mnl(x,fnl,schnl,dt,1);
ftlm = @(dx,x) Mtl(dx,x,ftl,fnl,schtl,schnl,dt,1); 

% Correctness Test for Lnl and Ltl
for i=1:16
alpha = 10^(1-i);
dx = alpha*randn(N,n+1);
x = randn(N,n+1);
xdx = x + dx;

Ldp = Loptl(dx,x,ftlm);
lpdp = Lopnl(xdx,fnlm);
lp   = Lopnl(x,fnlm);

alphadx(i)  = norm(lpdp - lp)/norm(Ldp);
alphadx2(i) = alphadx(i)-1;
alphadx3(i) = norm(lpdp - lp - Ldp);
end

semilogy(abs(alphadx)); figure
semilogy(abs(alphadx2)); figure 
semilogy(abs(alphadx3)); figure

% Correctness Test for Lnlinv and Ltlinv
for i=1:16
alpha = 10^(1-i);
dp = alpha*randn(N,n+1);
p = randn(N,n+1);
pdp = p + dp;

lpdp = Lopnlinv(pdp,fnlm);
lp   = Lopnlinv(p,fnlm);
Ldp  = Loptlinv(dp,lp,ftlm);

alphadxa(i)  = norm(lpdp - lp)/norm(Ldp);
alphadxa2(i) = alphadx(i)-1;
alphadxa3(i) = norm(lpdp - lp - Ldp);
end

semilogy(abs(alphadxa)); figure
semilogy(abs(alphadxa2)); figure 
semilogy(abs(alphadxa3))
