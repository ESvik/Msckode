%Richards Eq in 2D

%% Mesh
x_min=0; x_max=1; y_min=0; y_max=1;
h=1/8;
[x,y]=meshgrid(x_min:1/8:x_max,y_min:1/8:y_max);
X=reshape(x,[],1); Y=reshape(y,[],1);
DT=delaunayTriangulation(X,Y);
Coordinates=DT.Points;
Elements=DT.ConnectivityList;
Neumann=[];
Dirichlet=[find(Coordinates(:,1)==x_min);find(Coordinates(:,1)==x_max);find(Coordinates(:,2)==y_min);find(Coordinates(:,2)==y_max)];
DirichletValue=0;
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
    GaussValues(i,1) = 1/6*(Coefficients(i,1) + Coefficients(i,2)*1/6 + Coefficients(i,3)*1/6);
    GaussValues(i,2) = 1/6*(Coefficients(i,1) + Coefficients(i,2)*2/3 + Coefficients(i,3)*1/6);
    GaussValues(i,3) = 1/6*(Coefficients(i,1) + Coefficients(i,2)*1/6 + Coefficients(i,3)*2/3);
end

%% Time
t_0=7.9;
tau=0.1;
T=8;

%% Problem
uexact = @(x,y,t) x.*y.*(1-x).*(1-y).*t;
%Coefficients
k=0.01; L_b=10; steps=1; b_m=1;
pol=@(u) richardnonlin(u);
b= @(u) b_m.*u + (L_b-b_m)*(20/27)*pol(u);
%b = @(u) step(u,steps,L_b,b_m); 
f1= @(x,y,t) (b_m+(L_b-b_m)*(20/27).*(-5.4*uexact(x,y,t).^2+5.4*uexact(x,y,t))).*uexact(x,y,1)+2*k*t.*((1-y).*y+(1-x).*x);
f2= @(x,y,t) f1(x,y,t)*tau;
u_0 = @(x,y) 0;
g=0;
LL = @(g) 1./(2*(1-g));


%% Solver
L=L_b*0.5;
t=t_0+tau;
u = uexact(X,Y,t_0);
while t < T+tau
    f= @(x,y) f2(x,y,t);
    u_old=u;
    error = 1;
    while error > 10^(-8)
        u_prev=u;
        u=FEMParabolic2DP1Mscheme(Coordinates,Elements,L,k*tau,Dirichlet,DirichletValue,Neumann,g,f,zeros(2*NN,1),b(u_old)-b(u_prev)+L*u_prev);
        error=norm(u-u_prev,inf)+norm(u-u_prev,inf)/norm(u,inf)
    end
    t=t+tau;
end
uprecomp=u;
t=t_0+tau;
u = uexact(X,Y,t_0);
ErrorComparison = [];
while t < T+tau
    f= @(x,y) f2(x,y,t);
    u_old=u;
    error = 1; CompErrPrev=1;
    while error > 10^(-8)
        u_prev=u;
        u=FEMParabolic2D(Coordinates,Elements,L,k*tau,Dirichlet,DirichletValue,Neumann,g,f,zeros(2*NN,1),b(u_old)-b(u_prev)+L*u_prev);
        error=norm(u-u_prev,inf);
        CompErr=norm(u-uprecomp,inf);
        ErrorComparison = [ErrorComparison,CompErr/CompErrPrev];
        CompErrPrev=CompErr;
    end
    t=t+tau;
end


plot(ErrorComparison(2:end))
hold on
%end

% plot(Analysis(1,1:18)/L_b,Analysis(2,1:18)','k:')
% hold on
% plot(Analysis(1,19)/L_b,Analysis(2,19),'p','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',10)

% plot(Analysis(1,1:15)/L_b,Analysis(2,1:15)',':')
% hold on
%plot(Analysis(1,16)/L_b,Analysis(2,16),'p','MarkerEdgeColor','k','MarkerSize',10)


% legend('b_m=0.01','b_m=0.1','b_m=0.5','b_m=0.7','b_m=0.9')
%trisurf(Elements,X,Y,u)
% trisurf(Elements,X,Y,uexact(X,Y,T))