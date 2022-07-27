function [ ww,GG] = Forward_seismic_rock(m0,m_matrix,m_fluid,r,nt,uu,w0,ntheta,dtheta)

nt0=nt-1;
%%
D=DDchafen(nt);
DD1=zeros(nt0*2,nt*2);
DD1(1:nt0,1:nt)=D;
DD1(nt0+1:2*nt0,1+nt:2*nt)=D;

%%
SS=zeros(nt0*ntheta,nt0*2);
for itheta=1:ntheta  
    angle=(itheta-1)*dtheta+7.5;
    for i_r=1:nt0
        r_ratio=r(i_r,1);
        [AA]=Gassmann_model(angle,r_ratio,m0,m_matrix,m_fluid);
        a_por(i_r)=AA(1,1);
        a_cl(i_r)=AA(1,2);
    end
    ww=wave1(nt,w0,itheta); 
    ww=ww(1:nt0,1:nt0);
    iio=(itheta-1)*nt0;
        SS(iio+1:iio+nt0,1:nt0)=ww*diag(a_por);
        SS(iio+1:iio+nt0,1+nt0:2*nt0)=ww*diag(a_cl);
end
GG=SS*DD1*uu;

end

