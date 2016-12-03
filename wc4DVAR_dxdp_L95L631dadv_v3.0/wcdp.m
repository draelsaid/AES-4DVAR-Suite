% Weak-Constraint 4DVAR (dp/eta Formulation)

function [DJ,J] = wcdp(p,pb,y,fH,fHT,fD,fR,Mnl,Madj)
addpath(pwd,'L')

N = size(pb,1); % Number of variables determined by row size of p
s = size(pb,2); % Number of subwindows determined by col size of p

p = reshape(p,N,s);

pmpb = p - pb;                        % p-pb
Dinv_pmpb = fD(pmpb);                 % Dinv*(p - pb)
Jbq = sum(sum(pmpb.*Dinv_pmpb));      % Jb and Jq term combined

x = Lopnlinv(p,Mnl);                  % Linv(p)
h_Linv_p = fH(x);                     % H(Linv(p))
hlpmy = h_Linv_p - y;                 % Innovations
Rinv_hlpmy = fR(hlpmy);               % R Weighted Innovations

Jo = sum(sum(hlpmy.*Rinv_hlpmy));     % Jo term

J = 0.5*(Jbq + Jo);                   % Cost function

% Gradient Calc

LHRhlpmy = LoptlinvT(fHT(Rinv_hlpmy),x,Madj);

DJ = Dinv_pmpb + LHRhlpmy;            % Gradient
DJ = reshape(DJ,s*N,1);
end