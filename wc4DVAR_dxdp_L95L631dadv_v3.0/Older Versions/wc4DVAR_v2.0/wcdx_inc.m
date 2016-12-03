% Weak-Constraint 4DVAR (Incremental dx/x Formulation)

function [J,DJ,dx] = wcdx_inc(dx,x,b,d,fH,fHT,fD,fR,Mtl,Madj)
addpath(pwd,'L')

N = size(b,1); % Number of variables determined by row size of pb
s = size(b,2); % Number of subwindows determined by col size of pb

x = reshape(x,N,s); dx = reshape(dx,N,s);

Ldx = Loptl(dx,x,Mtl);                  % Ldx
Ldxb = Ldx-b;                           % Ldx - b
Dinv_Ldxb = fD(Ldxb);                   % Dinv*(Ldx-b)
Jbq = 0.5*sum(sum(Ldxb.*Dinv_Ldxb));    % Jb and Jq term combined

h_x = fH(dx);                           % Hdx
hxd = h_x - d;                          % Innovations
Rinv_hxd = fR(hxd);                     % R Weighted Innovations

Jo = 0.5*sum(sum(hxd.*Rinv_hxd));       % Jo term

J = Jbq + Jo;                           % Cost function

HRhxd = fHT(Rinv_hxd);                  % H^T Rinv (Hdx - d)
L_Dinv_Ldxb = LoptlT(Dinv_Ldxb,x,Madj); % L^T Dinv (Ldx - b)

DJ = L_Dinv_Ldxb + HRhxd;               % Cost function gradient 

DJ = reshape(DJ,s*N,1);
end
