%Richards Eq in 2D Nonlinear Diffusion

%% Mesh
x_min=0; x_max=1; y_min=0; y_max=1;
h=1/16;
[x,y]=meshgrid(x_min:h:x_max,y_min:h:y_max);
X=reshape(x,[],1); Y=reshape(y,[],1);
DT=delaunayTriangulation(X,Y);
Coordinates=DT.Points;
Elements=DT.ConnectivityList;
Neumann=[];
Dirichlet=zeros(4*(1/h+1),2);
Dirichlet(:,1)=[find(Coordinates(:,1)==x_min);find(Coordinates(:,1)==x_max);find(Coordinates(:,2)==y_min);find(Coordinates(:,2)==y_max)];
NN=length(Coordinates(:,1));

Coefficients= zeros(3,3);
MatOfCoord=[1,0,0;1,1,0;1,0,1];
for i =1:3
    bb=zeros(3,1);
    bb(i)=1;
    Coefficients(i,:)=MatOfCoord\bb;
end
GaussValues = zeros(3,3);
for i = 1:3
    GaussValues(i,1) = Coefficients(i,1) + Coefficients(i,2)*1/6 + Coefficients(i,3)*1/6;
    GaussValues(i,2) = Coefficients(i,1) + Coefficients(i,2)*2/3 + Coefficients(i,3)*1/6;
    GaussValues(i,3) = Coefficients(i,1) + Coefficients(i,2)*1/6 + Coefficients(i,3)*2/3;
end

%% Time
t_0=7.9;
tau=0.1;
T=t_0+tau;

%% Problem
uexact = @(x,y,t) x.*y.*(1-x).*(1-y).*t;
%Coefficients

kappa=@(u) u.^2+1;
kappaprime=@(u) 2*u;
L_b=10; b_m=1;
pol=@(u) richardnonlin(u);
b= @(u) b_m.*u + (L_b-b_m)*(20/27)*pol(u);
%b = @(u) step(u,steps,L_b,b_m); 
bprime = @(u) b_m+(L_b-b_m)*(20/27)*(-5.4*u.^2+5.4*u);
f1= @(x,y,t) (b_m+(L_b-b_m)*(20/27).*(-5.4*uexact(x,y,t).^2+5.4*uexact(x,y,t))).*uexact(x,y,1)+2*t.*((1-y).*y+(1-x).*x);
f2= @(x,y,t) f1(x,y,t)*tau;
g=0;

bprimemax=max(bprime(0:0.01:1));


%% Solver
KIRCHHOFF = 2;
SCHEME=1;
t=t_0+tau;
u_0 = uexact(X,Y,t_0);
f=@(x,y) f2(x,y,t);

a=10; beta=4*1*b_m; c=2*b_m; d=2*0.1; e=0.1;
gammaopt=max(0,(sqrt(beta*(a^2*beta+a*beta*(d-2*e)+a*d*(d+c)+(d-e)*(-beta*e+d*c)))+a*beta+beta*d-beta*e)/(beta*d));

if gammaopt>1-tau*(L_b)/2;
    gammaopt=0;
end
Analysis=zeros(100,2);
counter=1;
for gamma=[-0.3:0.05:7*0.05,gammaopt]

L=L_b/(2*(1-gamma)-tau*L_b);
[u,Au]=FEMParabolic2DP1Richards(Coordinates,Elements,L,0,Dirichlet,Neumann,g,f,b(u_0)-b(u_0),1,0,Coefficients,GaussValues,bprime,u_0,SCHEME,bprimemax/2,KIRCHHOFF,kappa,kappaprime,tau);
u_old = u_0;
error=norm(u-u_old,inf)+norm(u-u_old,inf)/norm(u,inf);
iterations=1;
while error > 10^(-8)
        u_prev=u;
        [u,Au]=FEMParabolic2DP1Richards(Coordinates,Elements,L,0,Dirichlet,Neumann,g,f,b(u_old)-b(u_prev),1,0,Coefficients,GaussValues,bprime,u_prev,SCHEME,bprimemax/2,KIRCHHOFF,kappa,kappaprime,tau);
        error=norm(u-u_prev,inf)+norm(u-u_prev,inf)/norm(u,inf)
        iterations=iterations+1;
end
Analysis(counter,1)=L/L_b; Analysis(counter,2)=iterations;
counter=counter+1;
end
hold on
plot(Analysis(1:14,1),Analysis(1:14,2),'b')
plot(Analysis(15,1),Analysis(15,2),'p','MarkerSize',10,'MarkerEdge','k','MarkerFace','k')
set(gca,'fontsize',25)


% end
% analysis(counter,2*index-1)=1.1; analysis(counter,2*index)=iterations;
% counter=counter+1;
% gammaopt=max(0,1+(L_b*b_m-sqrt(L_b^2*b_m^2+4*L_b*b_m^2*kappa*tau+4*L_b*b_m*kappa^2*tau^2))/(4*b_m*kappa*tau));
% 
% for gamma=[0:0.05:0.5,gammaopt]
%     kappa
%     gamma
% L=L_b/(2*(1-gamma));
% t=t_0+tau;
% u_0 = uexact(X,Y,t_0);
% f=@(x,y) f2(x,y,t);
% [u,Au]=FEMParabolic2DP1Richards(Coordinates,Elements,L,kappa*tau,Dirichlet,DirichletValue,Neumann,g,f,b(u_0)-b(u_0)+L*u_prev,1,0,Coefficients,GaussValues);
% u_old = u_0;
% error=norm(u-u_old,inf)+norm(u-u_old,inf)/norm(u,inf);
% iterations=1;
% while error > 10^(-8)
%         u_prev=u;
%         [u,Au]=FEMParabolic2DP1Richards(Coordinates,Elements,L,kappa*tau,Dirichlet,DirichletValue,Neumann,g,f,b(u_old)-b(u_prev)+L*u_prev,0,Au,Coefficients,GaussValues);
%         error=norm(u-u_prev,inf)+norm(u-u_prev,inf)/norm(u,inf);
%         iterations=iterations+1;
% end
% t=t+tau;
% 
% while t < T+tau
%     f= @(x,y) f2(x,y,t);
%     u_old=u;
%     error = 1;
%     while error > 10^(-8)
%         u_prev=u;
%         [u,Au]=FEMParabolic2DP1Richards(Coordinates,Elements,l,kappa*tau,Dirichlet,DirichletValue,Neumann,g,f,b(u_old)-b(u_prev)+L*u_prev,0,Au,Coefficients,GaussValues);
%         error=norm(u-u_prev,inf)+norm(u-u_prev,inf)/norm(u,inf);
%         iterations=iterations+1;
%     end
%     t=t+tau;
% end
% analysis(counter,2*index-1)=L/L_b; analysis(counter,2*index)=iterations;
% counter=counter+1;
% end
% 
% end
% plot(Analysis(1,1:18)/L_b,Analysis(2,1:18)','k:')
% hold on
% plot(Analysis(1,19)/L_b,Analysis(2,19),'p','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',10)
% plot(Analysis(1,1:15)/L_b,Analysis(2,1:15)',':')
% hold on
%plot(Analysis(1,16)/L_b,Analysis(2,16),'p','MarkerEdgeColor','k','MarkerSize',10)
% legend('b_m=0.01','b_m=0.1','b_m=0.5','b_m=0.7','b_m=0.9')
%trisurf(Elements,X,Y,u)