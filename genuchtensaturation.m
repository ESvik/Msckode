function s=genuchtensaturation(p,a,n)
    s=p;
    for i=1:length(p)
        if p(i)<=0
        s(i)=(1+(-a*p(i)).^n).^((-n+1)/n);
        end
        if p(i)>0 
        s(i)=1;
        end
    end
end