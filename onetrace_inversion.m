function [ m_est,Cm_est,cha ] = onetrace_inversion( miu,gg,Cm,Ce,seisd )
juzhen=inv(gg*Cm*gg'+Ce);
juzhen1=Cm*gg'*juzhen;
cha=juzhen1*(seisd-gg*miu);
lmiu=miu+cha;
Cm_est=Cm-Cm*gg'*juzhen*gg*Cm;
m_est(:,1)=lmiu;
end

