function [ rh1 ] = rh(c0,c,a,h,index1)

if index1==1
if h==0
    rh1=0;
else if h>a
        rh1=c0+c;
    else
        rh1=c0+c*((3*h)/(2*a)-(h*h*h)/(2*a*a*a));
    end
end
end

if index1==2
  if h==0
    rh1=0;
  else 
        a1=a/3;
        rh1=c0+c*(1-exp(-h/a1));
    end
end  

if index1==3
  if h==0
    rh1=0;
  else 
      a1=(a*a)/3;
        rh1=c0+c*(1-exp(-(h*h)/(a1*a1)));
    end
end  


if index1==4
     if h==0
    rh1=0;
  else 
      rh1=c0+c*(1-exp(-h/a)*cos((2*pi*h)/a));
    end
end  
end



