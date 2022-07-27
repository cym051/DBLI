function [ D ] = DDchafen( n )
D=diag(-ones(n-1,1),0)+diag(ones(n-2,1),1);
D=[D,zeros(n-1,1)];
D(end,end)=1;
end
