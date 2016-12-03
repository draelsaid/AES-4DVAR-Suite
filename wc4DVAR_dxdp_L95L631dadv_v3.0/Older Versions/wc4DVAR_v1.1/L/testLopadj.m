% Adjoint model test with numerical integration scheme
clear all
close all

addpath( pwd,'L')

N=40; dt=0.025; n=200;

fnl = @l95;
schnl = @rk4nl;

ftl = @l95tl;
schtl = @rk4tl;

fadj = @l95adj;
schadj = @rk4adj;

fLtl = @(dx,x) Lopl(dx,x,n,ftl,fnl,schtl,schnl,dt);
fLadj= @(dy,x) LoplT(dy,x,n,fadj,fnl,schadj,schnl,dt);


% Initial conditions 
x = randn(N,n+1); dx = randn(N,n+1); dy = randn(N,n+1);

%%%%%%%%%%%%%%%%%%%%%%%% Tests model with propogator (scheme) %%%%%%%%%%%%%%
Ldx = fLtl(dx,x);
Ltdy = fLadj(dy,x);

fwd_prd = sum(sum(dy.*Ldx)); disp(sprintf('%5.15f',fwd_prd))
adj_prd = sum(sum(Ltdy.*dx)); disp(sprintf('%5.15f',adj_prd))
Accuracy = fwd_prd - adj_prd
% Accuracy should be ~ 10^(-16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fLtl = @(dx,x) Loplinv(dx,x,n,ftl,fnl,schtl,schnl,dt);
fLadj= @(dy,x) LoplinvT(dy,x,n,fadj,fnl,schadj,schnl,dt);

%%%%%%%%%%%%%%%%%%%%%%%% Tests model with propogator (scheme) %%%%%%%%%%%%%%
Ldx = fLtl(dx,x);
Ltdy = fLadj(dy,x);

fwd_prd = sum(sum(dy.*Ldx)); disp(sprintf('%5.15f',fwd_prd))
adj_prd = sum(sum(Ltdy.*dx)); disp(sprintf('%5.15f',adj_prd))
Accuracy = fwd_prd - adj_prd
% Accuracy should be ~ 10^(-16)
