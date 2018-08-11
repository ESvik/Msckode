%Analytic mandel

function [p,ux,uy] = MandelAnalytic(F,B,vu,a,cf,mu,v,NumberOfSteps)
z=zeros(NumberOfSteps,1);
f=@(x) tan(x)-(1-v)/(vu-v)*x;
x0=[pi/4,pi/2-10000*eps];
options = optimset('TolFun',10^(-30),'TolX',10^-30);
for i=1:NumberOfSteps
[z(i),fval(i),info(i)] = fzero(f,x0,options);
x0=x0+pi;
end
%p=@(x,t) 2*F*B*(1+vu)/(3*a)*sum(sin(z)./(z-sin(z).*cos(z)).*(cos(z*x/a)-cos(z)).*exp(-z.^2*cf*t/a^2));
p=@(x,t)(( (2*F*B*(1+vu))/(3*a))*sum((sin(z)./(z-(sin(z).*cos(z)))).*(cos(z*x/a)-cos(z)).* exp(-(z.^2).*cf.*t./(a^2))));

ux=@(x,t) (F*v/(2*mu*a)-F*vu/(mu*a)*(sum(sin(z).*cos(z)./(z-sin(z).*cos(z)).*exp(-z.^2*cf*t/a^2)))).*x+F/mu*sum(cos(z)./(z-sin(z).*cos(z)).*sin(z*x/a).*exp(-z.^2*cf*t/a^2));
uy=@(y,t) (-F*(1-v)/(2*mu*a)+F*(1-vu)/(mu*a)*sum(sin(z).*cos(z)./(z-sin(z).*cos(z)).*exp(-z.^2*cf*t/a^2))).*y;
end