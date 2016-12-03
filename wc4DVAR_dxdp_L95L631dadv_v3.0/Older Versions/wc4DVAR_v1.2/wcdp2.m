% Weak-Constraint 4DVAR (dp/eta Formulation)

function [J,DJ] = wcdp(p,pb,y,H,s,n,fD,fR,fnl,fadj,schnl,schadj)

addpath( pwd,'L')

reshape(p,length(p)*n,1); reshape(pb,length(pb)*n,1);
pmpb = p-pb;
Jbq = pmpb' * Dinv * pmpb;

x = Lopnlinv(p,s,n,fnl,schnl,dt);
h_Linv_p = H*x;
reshape(y,length(y)*n,1); reshape(h_Linv_p,length(h_Linv_p)*n,1);
Jo = (y - h_Linv_p)' * Rinv * (y - h_Linv_p);

J = Jbq + Jo;

innov = H' * Rinv * (y - h_Linv_p);
LHRyhL = LoplinvT(innov,x,s,n,fadj,fnl,schadj,schnl,dt); 

DJ = Dinv * pmpb +  LHRyhL;



end
