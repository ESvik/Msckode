%%Newton solver

function x=Newton(x_0,f,fprime,tol,x_min,x_max)
x_old=x_0;
x_new=x_old-f(x_old)/fprime(x_0);
if x_new>x_max
    x_new=x_max;
end
if x_new<x_min
    x_new=x_min;
end
epsmax=(x_max-x_min)/10;
epsmin=epsmax;
while abs(x_new-x_old) > tol
    x_old=x_new;
    x_new=x_old-f(x_old)/fprime(x_old);
    if x_new>x_max
        x_new=x_max-epsmax;
        epsmax=epsmax+2*tol;
    end
    if x_new<x_min
        x_new=x_min+epsmin;
        epsmin=epsmin+2*tol;
    end
    x_new
end
x=x_new;
end
