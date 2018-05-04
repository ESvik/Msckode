% An increasing function in u with range from 0.1 to 1, with lipschitz
% constant 

function b=richardnonlin(u)
a=5.4;
b=u;
for i=[1:length(u)]
b(i)=-a/3*u(i)^3+a/2*u(i)^2+0.1;
if u(i)<=0
    b(i)=0.1;
end
if u(i)>=1
    b(i)=1;
end
end
end