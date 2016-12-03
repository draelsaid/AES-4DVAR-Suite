% Strong-Constraint 4DVAR
function [J,DJ] = sc4dvar(p,pb,y,fH,fHT,fD,fR,Mnl,Madj)
addpath(pwd,'L')

p = p(:,1); pb = pb(:,1);

N = size(pb,1); % Number of variables determined by row size of p
s = size(pb,2); % Number of subwindows determined by col size of p

p = reshape(p,N,s);

pmpb = p-pb;
Dinv_pmpb = fD(pmpb);                % Dinv*(p-pb)
Jb = 0.5*sum(sum(pmpb.*Dinv_pmpb));  % Jb and Jq term combined

x = Lopnlinv(p,Mnl);                 % Linv(p)
h_Linv_p = fH(x);                    % H(Linv(p))
yhlp = h_Linv_p - y;                 % Innovations
Rinv_yhlp = fR(yhlp);                % R Weighted Innovations

Jo = 0.5*sum(sum(yhlp.*Rinv_yhlp));  % Jo term

J = Jb + Jo;                         % Cost function

% Gradient Calc

LHRyhL = LoptlinvT(fHT(Rinv_yhlp),x,Madj);

DJ = Dinv_pmpb + LHRyhL;             % Gradient
DJ = reshape(DJ,s*N,1);
end
