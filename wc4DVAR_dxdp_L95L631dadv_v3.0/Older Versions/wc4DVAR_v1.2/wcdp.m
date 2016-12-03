% Weak-Constraint 4DVAR (dp/eta Formulation)

function [J,DJ] = wcdp(p,pb,y,n,fH,fD,fR,fnl,fadj,schnl,schadj,dt)
addpath(pwd,'L')
N = size(p,1); % Number of variables (or space) determined by row size of p
s = size(p,2); % Number of subwindows determined by col size of p

pmpb = p-pb;
Dinv_pmpb = fD(pmpb);
Jbq = sum(sum(pmpb.*Dinv_pmpb)); % Jb and Jq term combined

x = Lopnlinv(p,n,fnl,schnl,dt);  % Linv(p)
h_Linv_p = fH(x);                % H applied to Linv(p)
yhlp = y - h_Linv_p;             % Innovations
Rinv_yhlp = fR(yhlp);            % Weighted Innovations

Jo = sum(sum(yhlp.*Rinv_yhlp));  % Jo term

J = Jbq + Jo;

LHRyhL = LoplinvT(Rinv_yhlp,x,n,fadj,fnl,schadj,schnl,dt);

DJ = Dinv_pmpb +  LHRyhL;

end
