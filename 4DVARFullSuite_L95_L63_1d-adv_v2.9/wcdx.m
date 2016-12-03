% Weak-Constraint 4DVAR (dx/x Formulation)

function [DJ,J] = wcdx(x,pb,y,fH,fHT,fD,fR,Mnl,Madj)
addpath(pwd,'L')

N = size(pb,1); % Number of variables determined by row size of pb
s = size(pb,2); % Number of subwindows determined by col size of pb

x = reshape(x,N,s);

Lx = Lopnl(x,Mnl);                          % Lx
Lxmpb = Lx - pb;                            % Lx - pb
Dinv_Lxmpb = fD(Lxmpb);                     % Dinv*(Lx-pb)
Jbq = sum(sum(Lxmpb.*Dinv_Lxmpb));          % Jb and Jq term combined

h_x = fH(x);                                % h(x)
hxmy = h_x - y;                             % Innovations
Rinv_hxmy = fR(hxmy);                       % R Weighted Innovations

Jo = sum(sum(hxmy.*Rinv_hxmy));             % Jo term

J = 0.5*(Jbq + Jo);                         % Cost function

HRhxy = fHT(Rinv_hxmy);                     % H^T Rinv (h(x) - y)
L_Dinv_Lxmpb = LoptlT(Dinv_Lxmpb,x,Madj);   % L^T Dinv (L(x) - pb)
DJ = L_Dinv_Lxmpb + HRhxy;                  % Cost function gradient 

DJ = reshape(DJ,s*N,1);
end
