function [AA]=Gassmann_model(theta,r,m,m_matrix,m_fluid)
%%
%G=WSDm=WAFDm
%%
por=m(1,1);
cl=m(1,2);%

Kc=m_matrix(1,1);
Gc=m_matrix(1,2);
dc=m_matrix(1,3);
Kq=m_matrix(1,4);
Gq=m_matrix(1,5);
dq=m_matrix(1,6);
porc=0.38;
Kf=m_fluid(1,1);
df=m_fluid(1,2);


Km=(cl*Kc+(1-cl)*Kq+1/(cl/Kc+(1-cl)/Kq))/2;
Gm=(cl*Gc+(1-cl)*Gq+1/(cl/Gc+(1-cl)/Gq))/2;
dm=dc*cl+dq*(1-cl);
Ks=Km*(1-por/porc)+Kf*por/(porc^2+Km*Kf*porc);
Gs=Gm*(1-por/porc);

d=(1-por)*dm+por*df;
Ms=Ks+4*Gs/3;
vp=(Ms/d)^0.5;
vs=(Gs/d)^0.5;


d_por=df-dm;
Ks_por=(-Km+Kf/(porc+Km*Kf))/porc;
Gs_por=-Gm/porc;
A_por=((Ks_por+4*Gs_por/3)*d-Ms*d_por)/(d^2);
B_por=(Gs_por*d-Gs*d_por)/(d^2);


Km_cl=(Kc-Kq-(1/Kc-1/Kq)/(cl/Kc+(1-cl)/Kq)^2)/2;
Gm_cl=(Gc-Gq-(1/Gc-1/Gq)/(cl/Gc+(1-cl)/Gq)^2)/2;
dm_cl=dc-dq;
d_cl=(1-por)*dm_cl;
Ks_cl=Km_cl*(1-por/porc-(por*porc*Kf^2)/((porc^2+Km*Kf*porc)^2));
Gs_cl=Gm_cl*(1-por/porc);
A_cl=((Ks_cl+4*Gs_cl/3)*d-Ms*d_cl)/(d^2);
B_cl=(Gs_cl*d-Gs*d_cl)/(d^2);


vp_por=A_por/(2*vp);
vp_cl=A_cl/(2*vp);

vs_por=B_por/(2*vs);
vs_cl=B_cl/(2*vs);

angle=theta*pi/180;
A1=0.5*((sec(angle))^2);
A2=-4*(sin(angle))^2/(r^2);
A3=(1-4*(sin(angle))^2/(r^2))/2;

a_por=A1*vp_por/vp+A2*vs_por/vs+A3*d_por/d;
a_cl=A1*vp_cl/vp+A2*vs_cl/vs+A3*d_cl/d;

AA=[a_por a_cl];
end
