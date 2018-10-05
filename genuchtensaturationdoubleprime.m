function s=genuchtensaturationdoubleprime(p,a,n)
    s=p;
    for i=1:length(p)
        if p(i)<=0
        s(i)=a^2*(-n+1)*(1+(-a*p(i)).^n).^((-n+1)/n-1)+a^2*(-n+1)*((-n+1)/n-1)*p(i)*(1+(-a*p(i)).^n).^((-n+1)/n-2)*n*(-a*p(i))*(-a);
        end
        if p(i)>0 
        s(i)=0;
        end
    end
end