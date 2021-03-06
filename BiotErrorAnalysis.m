%% Biot Error Analysis

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
%Basis and its evaluation points for P1
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
%Basis and its evaluation points for P2
CoefficientsP2= zeros(3,6);
MatOfCoord=[1,0,0,0,0,0;1,1,0,1,0,0;1,0,1,0,1,0;1,1/2,0,1/4,0,0;1,0,1/2,0,1/4,0;1,1/2,1/2,1/4,1/4,1/4];   
CoefficientsP2(1,:)=MatOfCoord\[1;0;0;1/2;1/2;0];
CoefficientsP2(2,:)=MatOfCoord\[0;1;0;1/2;0;1/2];
CoefficientsP2(3,:)=MatOfCoord\[0;0;1;0;1/2;1/2];    
CoefficientsP2=[CoefficientsP2(1,:);CoefficientsP2(1,:);CoefficientsP2(2,:);CoefficientsP2(2,:);CoefficientsP2(3,:);CoefficientsP2(3,:)];
GaussValuesP2 = zeros(6,3);
for i = 1:6
    GaussValuesP2(i,1) = (CoefficientsP2(i,1) + CoefficientsP2(i,2)*1/6 + CoefficientsP2(i,3)*1/6+CoefficientsP2(i,4)*1/6^2+CoefficientsP2(i,5)*1/6^2+CoefficientsP2(i,6)*1/6^2);
    GaussValuesP2(i,2) = (CoefficientsP2(i,1) + CoefficientsP2(i,2)*2/3 + CoefficientsP2(i,3)*1/6+CoefficientsP2(i,4)*(2/3)^2+CoefficientsP2(i,5)*1/6^2+CoefficientsP2(i,6)*1/6*2/3);
    GaussValuesP2(i,3) = (CoefficientsP2(i,1) + CoefficientsP2(i,2)*1/6 + CoefficientsP2(i,3)*2/3+CoefficientsP2(i,4)*1/6^2+CoefficientsP2(i,5)*(2/3)^2+CoefficientsP2(i,6)*1/6*2/3);
end
Dirichletp=[find(Coordinates(:,1)==x_min);find(Coordinates(:,1)==x_max);find(Coordinates(:,2)==y_min);find(Coordinates(:,2)==y_max)];
Dirichletu = zeros(8*(1/h+1),2);
Dirichletu(:,1)=[2*find(Coordinates(:,1)==x_min);2*find(Coordinates(:,1)==x_min)-ones(length(find(Coordinates(:,1)==x_min)),1);2*find(Coordinates(:,1)==x_max);2*find(Coordinates(:,1)==x_max)-ones(length(find(Coordinates(:,1)==x_max)),1);2*find(Coordinates(:,2)==y_min);2*find(Coordinates(:,2)==y_min)-ones(length(find(Coordinates(:,2)==y_min)),1);2*find(Coordinates(:,2)==y_max);2*find(Coordinates(:,2)==y_max)-ones(length(find(Coordinates(:,2)==y_max)),1)];
Dirichletu(1:2*(1/h+1),2)=0;
Dirichletu(2*(1/h+1)+1:4*(1/h+1),2)=0;
Dirichletu(4*(1/h+1)+1:6*(1/h+1),2)=0;
Dirichletu(6*(1/h+1)+1:8*(1/h+1),2)=0;

DirichletValue=0;

%% Time
t_0=0;
tau=10^(-1);
T=10^(-1);

%% Problem
pressurescale=10^13;
uexact = @(x,y,t) [t.*x.*y.*(x-1).*(y-1),t.*x.*y.*(x-1).*(y-1)];
pexact = @(x,y,t) pressurescale*t.*x.*y.*(x-1).*(y-1);
lambda = 27.778*10^(9); mu=41.667*10^(9); M=100*10^10; alpha=1; kappa=10^(-12);
%lambda = 1; mu=1; M=1; alpha=1; kappa=10^(0);
u_0=zeros(2*NN,1);
u0=uexact(X,Y,t_0);
u_0(1:2:2*NN-1)=u0(:,1);
u_0(2:2:2*NN)=u0(:,2);
p_0=pexact(X,Y,t_0);
%f_1=@(x,y,t) [0;0];
%f_2=@(x,y,t) 0;
f_1=@(x,y,t) [(-2*mu-lambda)*2*t*y.*(y-1)+(-mu-lambda)*(2*x-1).*(2*y-1).*t-mu*2*t*x.*(x-1)+(alpha*t*y.*(y-1).*(2*x-1))*pressurescale; -mu*2*t*y.*(y-1)+(-mu-lambda)*t*(2*x-1).*(2*y-1)+(-2*mu-lambda)*2*t*x.*(x-1)+(alpha*t*x.*(x-1).*(2*y-1))*pressurescale];
f_2=@(x,y,t) (1/M*x.*y.*(x-1).*(y-1)-kappa*t*2*(x.*(x-1)+y.*(y-1)))*pressurescale+alpha*(y.*(y-1).*(2*x-1)+x.*(x-1).*(2*y-1));

%% Mathematical optima
A_delta=(2/M+2*tau*kappa+2*alpha^2/(2*mu+lambda));
B_delta=(alpha^2/(2*mu+lambda));
delta_opt=min(A_delta/(2*B_delta),2);
%% Stabilization parameter
delta = 2;
L=alpha^2/((mu+lambda)*delta);
scale=1;
%% Precomputing
t=t_0+tau;


f_10=@(x,y) f_1(x,y,t);
f_20=@(x,y) tau*f_2(x,y,t);
[p,Ap] = FEMParabolic2D(Coordinates,Elements,1/M+L,kappa*tau,Dirichletp,DirichletValue,Neumann,g,f_20,zeros(2*NN,1),1/M*p_0+L*p_0,1,0,Coefficients,GaussValues);
[u,Au] = FEMLinStrain2DvectorialP2(Coordinates,Elements,2*mu,lambda,Dirichletu,Neumann,g,f_10,alpha*p,zeros(2*NN,1),1,0,CoefficientsP2,GaussValuesP2);
u_old = u_0;
p_old = p_0;
errorp=norm(p-p_old,inf)/norm(p,inf)
erroru=norm(u-u_old,inf)/norm(u,inf)
    while errorp>10^(-12) || erroru>10^(-12)
        u_prev=u;
        p_prev=p;
        [p,~] = FEMParabolic2D(Coordinates,Elements,1/M+L,kappa*tau,Dirichletp,DirichletValue,Neumann,g,f_20,-alpha*(u_prev-u_old),1/M*p_old+L*p_prev,0,Ap,Coefficients,GaussValues);
        [u,~] = FEMLinStrain2DvectorialP2(Coordinates,Elements,2*mu,lambda,Dirichletu,Neumann,g,f_10,alpha*p,zeros(2*NN,1),0,Au,CoefficientsP2,GaussValuesP2);
        errorp=norm(p-p_prev,inf)/norm(p,inf)
        erroru=norm(u-u_prev,inf)/norm(u,inf)
        %error = sqrt(errorp^2+erroru^2)
    end
t=t+tau



while t<T
    u_old = u;
    p_old = p;
    f_1t = @(x,y) f_1(x,y,t);
    f_2t = @(x,y) tau*f_2(x,y,t)*scale;
    errorp = 1; erroru=1;
    while errorp>10^(-12) || erroru>10^(-12)
        u_prev=u;
        p_prev=p;
        [p,~] = FEMParabolic2D(Coordinates,Elements,1/M+L,kappa*tau,Dirichletp,DirichletValue,Neumann,g,f_20,-alpha*(u_prev-u_old),1/M*p_old+L*p_prev,0,Ap,Coefficients,GaussValues);
        [u,~] = FEMLinStrain2DvectorialP2(Coordinates,Elements,2*mu,lambda,Dirichletu,Neumann,g,f_10,alpha*p,zeros(2*NN,1),0,Au,CoefficientsP2,GaussValuesP2);
        errorp=norm(p-p_prev,inf)/norm(p,inf);
        erroru=norm(u-u_prev,inf)/norm(u,inf);
        %error = sqrt(errorp^2+erroru^2)
    end
    t=t+tau
end
uPreComp=u;
pPreComp=p;

%% Error analysis
% ErrorComparisonP=[];
% ErrorComparisonU=[];
% ErrCompPrevP=1;
% ErrCompPrevU=1;
% t=t_0+tau;
% u=u_0;
% p=p_0;
% while t<T+tau
%     u_old = u;
%     p_old = p;
%     f_1t = @(x,y) f_1(x,y,t);
%     f_2t = @(x,y) tau*f_2(x,y,t);
%     errorp = 1; erroru=1;
%     while errorp>10^(-12) | erroru>10^(-12)
%         u_prev=u;
%         p_prev=p;
%         p=FEMParabolic2D(Coordinates,Elements,1/M+L,kappa*tau,Dirichletp,DirichletValue,Neumann,g,f_2t,-alpha*(u_prev-u_old),1/M*p_old+L*p_prev);
%         u=FEMLinStrain2DvectorialP2(Coordinates,Elements,2*mu,lambda,Dirichletu,DirichletValue,Neumann,g,f_1t,alpha*p,zeros(2*NN,1));
%         errorp=norm(p-p_prev,inf)/norm(p,inf); erroru=norm(u-u_prev,inf)/norm(u,inf);
%         %error = sqrt(errorp^2+erroru^2)
%         ErrCompP=norm(p-pPreComp,inf); ErrCompU=norm(u-uPreComp,inf);
%         ErrorComparisonP=[ErrorComparisonP,ErrCompP/ErrCompPrevP];
%         ErrorComparisonU=[ErrorComparisonU,ErrCompU/ErrCompPrevU];
%         ErrCompPrevP=ErrCompP;
%         ErrCompPrevU=ErrCompU;
%     end
%     t=t+tau;
% end

% subplot(2,1,1)
% plot(ErrorComparisonU(2:end))
% hold on
% subplot(2,1,2)
% plot(ErrorComparisonP(2:end))
% hold on

%% Plots
subplot(3,1,1)
trisurf(Elements,X,Y,uPreComp(1:2:2*NN-1))
subplot(3,1,2)
trisurf(Elements,X,Y,uPreComp(2:2:2*NN))
subplot(3,1,3)
trisurf(Elements,X,Y,pPreComp)