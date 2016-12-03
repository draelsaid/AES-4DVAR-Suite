% Weak-Constraint 4DVAR (Incremental dp/p Formulation)

function [DJ,J,dp] = wcdp_inc(dp,p,b,d,fH,fHT,fD,fR,Mtl,Madj)
addpath(pwd,'L')

N = size(b,1); % Number of variables determined by row size of pb
s = size(b,2); % Number of subwindows determined by col size of pb

p = reshape(p,N,s); dp = reshape(dp,N,s);

dpmb = dp-b;
Dinv_dpmb = fD(dpmb);                % Dinv*(dp-b)
Jbq = 0.5*sum(sum(dpmb.*Dinv_dpmb)); % Jb and Jq term combined

dp = Loptlinv(dp,p,Mtl);             % Linv(dp)
h_Linv_dp = fH(dp);                  % H(Linv(p))
yhldp = h_Linv_dp - d;               % Innovations
Rinv_yhldp = fR(yhldp);               % R Weighted Innovations

Jo = 0.5*sum(sum(yhldp.*Rinv_yhldp)); % Jo term

J = Jbq + Jo;                        % Cost function

LHRyhL = LoptlinvT(fHT(Rinv_yhldp),p,Madj);

DJ = Dinv_dpmb + LHRyhL;

DJ = reshape(DJ,s*N,1);
end
