clc;clear;

load ini_rock_resample;
por=ini_por(:,1:340);
cl=ini_cl(:,1:340);
vp=ini_vp(:,1:340);
vs=ini_vs(:,1:340);
r=vp./vs;

N=length(por(1,:));
M=length(por(:,1));
ntd=M-1;
por0=mean(mean(por));
cl0=mean(mean(cl)); 
m0=[por0,cl0];

% rock physics parameters
Kc=21;
Gc=10;
dc=2.1;
Kq=36;
Gq=44;
dq=2.65;
Kf=2.3;
df=1; 

m_matrix=[Kc Gc dc Kq Gq dq];
m_fluid=[Kf df];

load seisd_jiami;
seis7_5=seis(1:M,:);
seis17_5=seis(M+1:M*2,:);
seis27_5=seis(1+M*2:M*3,:);
for i=1:N
    seis7_5(:,i)=seis7_5(:,i)/max(abs(seis7_5(:,i)));
    seis17_5(:,i)=seis17_5(:,i)/max(abs(seis17_5(:,i)));
    seis27_5(:,i)=seis27_5(:,i)/max(abs(seis27_5(:,i)));   
end

%Wave
load sub_wave;
sub_wave=sub_wave/(max(abs(sub_wave)));
w0(:,1)=sub_wave;
w0(:,2)=sub_wave;
w0(:,3)=sub_wave;
ntheta=3;dtheta=10;

for i=1:N
    poryh(:,i)=smooth(por(1:M,i),20); 
    clyh(:,i)=smooth(cl(1:M,i),20);  
end
for j=1:M
    poryh(j,:)=smooth(poryh(j,:),20); 
    clyh(j,:)=smooth(clyh(j,:),20);  
end

%Geostatistical parameters
cpor0_hor=0;ccl0_hor=0;
cpor0_ver=0.000;ccl0_ver=0.000;

cpora_hor=0.01;ccla_hor=0.035;
cpora_ver=0.01;ccla_ver=0.035;

apor_ver=20;acl_ver=20;
apor_hor=20;acl_hor=20;

for itrace=1:N
   disp(itrace) 
   seis_d=zeros(ntd*ntheta,1);
   cc=zeros(M,2);
   
   seis_d(1:ntd,1)=seis7_5(1:ntd,itrace);
   seis_d(ntd+1:2*ntd,1)=seis17_5(1:ntd,itrace);
   seis_d(2*ntd+1:3*ntd,1)=seis27_5(1:ntd,itrace);
   cc(1:M,1)=poryh(:,itrace);
   cc(1:M,2)=clyh(:,itrace);
   rr(1:M,1)=r(:,itrace);
   
   %Seismic forward model
   [aaww,gg]=Forward_seismic_rock(m0,m_matrix,m_fluid,r,M,uu,w0,ntheta,dtheta);
   %% Decorrelation
   covc=cov(cc);
   [u,S,D] = svd(covc);
   v=u';
   uni=v;
   uu=zeros(M*2);
   for i=1:2
       for j=1:2
           uuu=ones(M,1)*u(i,j);
           uu((i-1)*M+1:i*M,(j-1)*M+1:j*M)=diag(uuu);
       end
   end
   
   porc=cc(:,1);
   clc=cc(:,2);

   cc1=zeros(M,2);
   cc1(1:M,1)=porc;
   cc1(1:M,2)=clc;
   
   czhuanzhi=uni*cc1';
   cc2=czhuanzhi';
   cc1=cc2;
   porc=cc1(:,1);
   clc=cc1(:,2);
   
   miu=zeros(M*2,1);   
   miu(1:ntn,1)=porc;
   miu(1+ntn:2*ntn,1)=clc;
      
   Cm=zeros(2*M,2*M);
   [cpor]=multrace_Cm(M,cpor0_ver,cpora_ver,apor_ver,cpor0_hor,cpora_hor,apor_hor);
   [ccl]=multrace_Cm(M,ccl0_ver,ccla_ver,acl_ver,ccl0_hor,ccla_hor,acl_hor);   
   Cm(1:M,1:M)=cpor;
   Cm(M+1:2*M,M+1:2*M)=ccl;
   
   ee=0.01*ones(size(seis_d));
   Ce=diag(ee);
   %% posterior mean
   [m_est,Cm_est,cha]=onetrace_inversion(miu,gg,Cm,Ce,seis_d);
   core_m_est=uu*m_est;
   core_Cm_est=uu*Cm_est*uu';
   
   i_por=core_m_est(1:M);
   i_cl=core_m_est(M+1:2*M);
   
   i_porc=core_Cm_est(1:M,1:M);
   i_clc=core_Cm_est(M+1:2*M,M+1:2*M);
   
   inv_por(:,itrace)=i_por;
   inv_cl(:,itrace)=i_cl;
end
