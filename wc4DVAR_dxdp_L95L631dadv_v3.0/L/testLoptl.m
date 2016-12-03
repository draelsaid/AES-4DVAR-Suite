% Testing the tangent-linear model with the numerical integration scheme
clear all
close all

addpath( pwd,'L')

N=50; dt=0.02; n=50; 

fnl = @adv1d; schnl = @lpde;
ftl = @adv1d; schtl = @lpdetl;

fnlm = @(x) Mnl(x,fnl,schnl,dt,1);
ftlm = @(dx,x) Mtl(dx,x,ftl,fnl,schtl,schnl,dt,1); 

% Correctness Test for Lnl and Ltl
for i=17:-1:1
    alpha(i) = 10^(1-i);
end
alpha(18) = 0;

for i=1:18
dx = alpha(i)*randn(N,n+1);
x = randn(N,n+1);
xdx = x + dx;

Ldx = Loptl(dx,x,ftlm);
lxdx = Lopnl(xdx,fnlm);
lx   = Lopnl(x,fnlm);

% ldx = Lopnl(dx,fnlm); 
% Lxdx = Loptl(xdx,x,ftlm);
% 
% test(i) = norm(ldx) / norm(Ldx);
% test2(i) = norm(lxdx) / norm(Lxdx);

alphadx(i)  = norm(lxdx - lx)/norm(Ldx);
alphadx2(i) = 1-alphadx(i);
alphadx3(i) = norm(lxdx - lx - Ldx);
end

loglog(alpha,abs(alphadx)); figure
loglog(alpha,abs(alphadx2)); figure 
loglog(alpha,abs(alphadx3)); figure

% Correctness Test for Lnlinv and Ltlinv
for i=1:18
dp = alpha(i)*randn(N,n+1);
p = randn(N,n+1);
pdp = p + dp;

lpdp = Lopnlinv(pdp,fnlm);
lp   = Lopnlinv(p,fnlm);
Ldp  = Loptlinv(dp,lp,ftlm);

alphadxa(i)  = norm(lpdp - lp)/norm(Ldp);
alphadxa2(i) = 1-alphadxa(i);
alphadxa3(i) = norm(lpdp - lp - Ldp);
end

loglog(alpha,abs(alphadxa)); figure
loglog(alpha,abs(alphadxa2)); figure
loglog(alpha,abs(alphadxa3))
