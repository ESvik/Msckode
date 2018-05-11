%% Biot Analysis
%% MESH
x_min=0; x_max=1; y_min=0; y_max=1;
h=0.1;
[x,y]=meshgrid(x_min:h:x_max,y_min:h:y_max);
X=reshape(x,[],1); Y=reshape(y,[],1);
DT=delaunayTriangulation(X,Y);
Coordinates=DT.Points;
Elements=DT.ConnectivityList;
Neumann=[];
g=0;
NN=length(X);
Dirichletp=[find(Coordinates(:,1)==x_min);find(Coordinates(:,1)==x_max);find(Coordinates(:,2)==y_min);find(Coordinates(:,2)==y_max)];
Dirichletu = zeros(8*(1/h+1),2);
Dirichletu(:,1)=[2*find(Coordinates(:,1)==x_min);2*find(Coordinates(:,1)==x_min)-ones(length(find(Coordinates(:,1)==x_min)),1);2*find(Coordinates(:,1)==x_max);2*find(Coordinates(:,1)==x_max)-ones(length(find(Coordinates(:,1)==x_max)),1);2*find(Coordinates(:,2)==y_min);2*find(Coordinates(:,2)==y_min)-ones(length(find(Coordinates(:,2)==y_min)),1);2*find(Coordinates(:,2)==y_max);2*find(Coordinates(:,2)==y_max)-ones(length(find(Coordinates(:,2)==y_max)),1)];


DirichletValue=0;

%% Time
t_0=0.25;
tau=0.01;
T=0.26;


%% Problem
uexact = @(x,y,t) [t.*x.*y.*(x-1).*(y-1),t.*x.*y.*(x-1).*(y-1)];
vexact = @(x,y,t) t.*x.*y.*(x-1).*(y-1);
lambda = 27.778*10^(9); mu=41.667*10^(9); M=100*10^10; alpha=1; kappa=10^(-12);
u_0=zeros(2*NN,1);
u0=uexact(X,Y,t_0);
u_0(1:2:2*NN-1)=u0(:,1);
u_0(2:2:2*NN)=u0(:,2);
p_0=vexact(X,Y,t_0);
%f_1=@(x,y,t) [0;0];
%f_2=@(x,y,t) 0;
pressurescale=10^10;
f_1=@(x,y,t) [(-2*mu-lambda)*2*t*y.*(y-1)+(-mu-lambda)*(2*x-1).*(2*y-1).*t-mu*2*t*x.*(x-1)+(alpha*t*y.*(y-1).*(2*x-1))*pressurescale; -mu*2*t*y.*(y-1)+(-mu-lambda)*t*(2*x-1).*(2*y-1)+(-2*mu-lambda)*2*t*x.*(x-1)+(alpha*t*x.*(x-1).*(2*y-1))*pressurescale];
f_2=@(x,y,t) (1/M*x.*y.*(x-1).*(y-1)-kappa*t*2*(x.*(x-1)+y.*(y-1)))*pressurescale+alpha*(y.*(y-1).*(2*x-1)+x.*(x-1).*(2*y-1));



%% Solver
Analysis=zeros(12,2);
A_delta=(2/M+2*tau*kappa+2*alpha^2/(2*mu+lambda));
B_delta=(alpha^2/(2*mu+lambda));
delta_opt=min(A_delta/(2*B_delta),2)
counter=1;
%for delta = [1.5:0.1:3,delta_opt]
    delta=4
    iterations=0;
    L=alpha^2/((mu+lambda)*delta);
    t=t_0+tau;
    u=u_0;
    p=p_0;
    while t<T+tau
        Dirichletu(1:2*(1/h+1),2)=0;
        Dirichletu(2*(1/h+1)+1:4*(1/h+1),2)=0;
        Dirichletu(4*(1/h+1)+1:6*(1/h+1),2)=0;
        Dirichletu(6*(1/h+1)+1:8*(1/h+1),2)=0;
        u_old = u;
        p_old = p;
        f_1t = @(x,y) f_1(x,y,t);
        f_2t = @(x,y) tau*f_2(x,y,t);
        erroru = 1; errorp = 1;
        while errorp>10^(-8) | erroru>10^(-8)
            u_prev=u;
            p_prev=p;
            p=FEMParabolic2D(Coordinates,Elements,1/M+L,kappa*tau,Dirichletp,DirichletValue,Neumann,g,f_2t,-alpha*(u_prev-u_old),1/M*p_old+L*p_prev);
            u=FEMLinStrain2DvectorialP2(Coordinates,Elements,2*mu,lambda,Dirichletu,DirichletValue,Neumann,g,f_1t,alpha*p,zeros(2*NN,1));
            errorp=norm(p-p_prev,inf)/norm(p,inf)
            erroru=norm(u-u_prev,inf)/norm(u,inf)
            iterations=iterations+1;
        end
        t=t+tau;
    end
    iterations
    Analysis(counter,:)=[delta,iterations];
    counter=counter+1;

%end
plot(Analysis(1:16,1),Analysis(1:16,2))
hold on
%plot(Analysis(17,1),Analysis(17,2),'p','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',10)

%subplot(3,1,1)
%trisurf(Elements,X,Y,u(1:2:2*NN-1))
%subplot(3,1,2)
%trisurf(Elements,X,Y,u(2:2:2*NN))
%subplot(3,1,3)
%trisurf(Elements,X,Y,p)