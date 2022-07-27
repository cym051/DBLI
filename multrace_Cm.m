function [cm ] = multrace_Cm( nt,c0_ver,ca_ver,a_ver,c0_hor,ca_hor,a_hor)
index1=2;
cm0=zeros(nt);
for i=1:nt
    for j=1:nt
        if j>i
            hh=(j-i);
        else
            hh=(i-j);
        end
        sv=rh(c0_ver,ca_ver,a_ver,hh,index1);
        cm0(i,j)=ca_ver-sv;
    end
end
cm=zeros(nt);
trace_location1=[1 2];
for i=1:1
    for j=1:1
        distance_col=abs(trace_location1(i)-trace_location1(j));
        Variogram_col=rh(c0_hor,ca_hor,a_hor,distance_col,index1);
        ca_col=ca_hor-Variogram_col;
        cm_col_ratio=ca_col/ca_hor;
        cm((i-1)*nt+1:i*nt,(j-1)*nt+1:j*nt)=cm_col_ratio*cm0;
    end
end
end

