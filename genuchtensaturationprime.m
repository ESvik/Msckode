function s=genuchtensaturationprime(p,a,n)
    s=p;
    for i=1:length(p)
        if p(i)<=0
        s(i)=a^2*(-n+1)*p(i)*(1+(-a*p(i)).^n).^((-n+1)/n-1);
        end
        if p(i)>0 
        s(i)=0;
        end
    end
end