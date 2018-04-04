%generates a step function between 0 and 1 with m steps, and maximal rate of change
% a>=1, and minimal rate of change b<=1 where a=1 or b=1 corresponds to
% fully linear case and no steps

function u = step(x,m,a,b)
m=m*2;
u=zeros(1,length(x));
for i=1:length(x)
if x(i)<0
    u(i)=0;
end
if x(i)>1
    u(i)=1;
end
for n=0:2/m:1
    if x(i)>=n
        if x(i)<=n+2/(m*(a-b))*(a-1)
            u(i)=n+b*(x(i)-n);
        end
    end
    if x(i)>=n+2/(m*(a-b))*(a-1)
        if x(i)<=n+2/m
            u(i)=a*(x(i)-n-2/m)+n+2/m;
        end
    end
end
end
end
    