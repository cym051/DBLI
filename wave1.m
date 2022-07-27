function [ ww1 ] = wave1(nt,www,i_ang)
wl=www(:,i_ang);
nwt=length(wl);
nt0=nwt/2;
w=zeros(nwt,nt);
for ii=1:nt
    w(:,ii)=wl;
end
ww1=zeros(nt);
for it1=1:nt
    for it2=1:nt
        tt0=it2-it1+nt0+1;
        if tt0<=nwt&&tt0>=1
            ww1(it1,it2)=w(tt0,it2);
        end
    end
end
end

