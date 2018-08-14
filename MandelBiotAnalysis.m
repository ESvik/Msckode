%% Biot convergence optimization analysis on Mandel's problem
%% MESH
x_min=0; x_max=100; y_min=0; y_max=10;
dx=12.5; dy=2.5;
Nx=(x_max-x_min)/dx+1; Ny=(y_max-y_min)/dy+1;
[x,y]=meshgrid(x_min:dx:x_max,y_min:dy:y_max);
X=reshape(x,[],1); Y=reshape(y,[],1);
DT=delaunayTriangulation(X,Y);
Coordinates=DT.Points;
Elements=DT.ConnectivityList;
Neumann=[];
g=0;
NN=length(X);
[ElementsP2,CoordinatesP2]=ElementsPlusEdgesMandel(Elements,Coordinates,dx,dy,x_min,x_max,y_min,y_max);
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
CoefficientsP2= zeros(6,6);
MatOfCoord=[1,0,0,0,0,0;1,1,0,1,0,0;1,0,1,0,1,0;1,1/2,0,1/4,0,0;1,0,1/2,0,1/4,0;1,1/2,1/2,1/4,1/4,1/4];   
CoefficientsP2(1,:)=MatOfCoord\[1;0;0;0;0;0];
CoefficientsP2(2,:)=MatOfCoord\[0;1;0;0;0;0];
CoefficientsP2(3,:)=MatOfCoord\[0;0;1;0;0;0];
CoefficientsP2(4,:)=MatOfCoord\[0;0;0;1;0;0];
CoefficientsP2(5,:)=MatOfCoord\[0;0;0;0;1;0];
CoefficientsP2(6,:)=MatOfCoord\[0;0;0;0;0;1];
CoefficientsP2=[CoefficientsP2(1,:);CoefficientsP2(1,:);CoefficientsP2(2,:);CoefficientsP2(2,:);CoefficientsP2(3,:);CoefficientsP2(3,:);CoefficientsP2(4,:);CoefficientsP2(4,:);CoefficientsP2(5,:);CoefficientsP2(5,:);CoefficientsP2(6,:);CoefficientsP2(6,:)];
GaussValuesP2 = zeros(12,3);
for i = 1:12
    GaussValuesP2(i,1) = (CoefficientsP2(i,1) + CoefficientsP2(i,2)*1/6 + CoefficientsP2(i,3)*1/6+CoefficientsP2(i,4)*1/6^2+CoefficientsP2(i,5)*1/6^2+CoefficientsP2(i,6)*1/6^2);
    GaussValuesP2(i,2) = (CoefficientsP2(i,1) + CoefficientsP2(i,2)*2/3 + CoefficientsP2(i,3)*1/6+CoefficientsP2(i,4)*(2/3)^2+CoefficientsP2(i,5)*1/6^2+CoefficientsP2(i,6)*1/6*2/3);
    GaussValuesP2(i,3) = (CoefficientsP2(i,1) + CoefficientsP2(i,2)*1/6 + CoefficientsP2(i,3)*2/3+CoefficientsP2(i,4)*1/6^2+CoefficientsP2(i,5)*(2/3)^2+CoefficientsP2(i,6)*1/6*2/3);
end
%Boundary for flow
Dirichletp=find(Coordinates(:,1)==x_max);
DirichletValue = 0;
% Neumannp = zeros(3*(1/h+1),2);
% Neumannp(:,1)=[find(Coordinates(:,1)==x_min);find(Coordinates(:,2)==y_min);find(Coordinates(:,2)==y_max)];
% Neumannp(:,2)=0;
%Boundary for mechanics
NxP2=length(find(CoordinatesP2(:,2)==y_max)); NyP2=length(find(CoordinatesP2(:,1)==x_min));
Dirichletu = zeros(2*NxP2+NyP2,2);
Dirichletu(:,1)=[2*find(CoordinatesP2(:,2)==y_max);2*find(CoordinatesP2(:,1)==x_min)-ones(length(find(CoordinatesP2(:,1)==x_min)),1);2*find(CoordinatesP2(:,2)==y_min)];
% Dirichletu(1:(10/dy+1),2)=0;
% Dirichletu(10/dy+2:10/dy+1+(100/dx+1),2)=0;



%% Problem
a = 100;
b = 10; 
E =  5.94*10^9;% [Pa]
v = 0.2;% 
alpha=1;% 
B=0.833333; %Skempton coefficient 
vu= 0.44; %Undrained Poison's ratio
M=1.65*10^10;% Biot Modulus [Pa]
F=6*10^8;
lambda=(E*v)/((1-2*v)*(1+v));
mu= E/(2+2*v);
K=lambda+2/3*mu;
Ku= K + alpha^2*M;
visc=0.01;% Viscocity

Beta=((1-2*v)*(1+v))/(2*E*(1-v));
kdr=1/Beta;
%mu=1;
%c=0.465; %[m2/s]
%k=100;
%c=(2*kappa*(B^2)*mu*(1-v)*(1+vu)^2) /(9*visc*(1-vu)*(vu-v));
PresNormalizer=2*F*B*(1+vu)/(3*a);
%lambda = 1; mu=1; M=1; alpha=1; kappa=10^(0);
u_0=zeros(2*(2*Nx-1)*(2*Ny-1),1);
p_0=zeros(NN,1);
f_1=@(x,y,t) [0;0];
f_2=@(x,y,t) 0;
Analysis=zeros(17,12);
kappavector=[10^(-14),10^(-13),10^(-12),10^(-11),10^(-10)];


%for index=1:6
   kappa = kappavector(5);
cf = M*kappa*(K+4/3*mu)/(Ku+4/3*mu);
%cf = (2*kappa*(B^2)*mu*(1-v)*(1+vu)^2) /(9*visc*(1-vu)*(vu-v));

%% Analytical solution

[pAnalytic,uxAnalytic,uyAnalytic]=MandelAnalytic(F,B,vu,100,cf,mu,v,300);


%% Time
t_0=0;
tau=10^(4);
T=5*tau;

%% Mathematical optima
Kdr=1.7*mu+lambda;
beta=Kdr;
A_delta=(2/M+2*tau*kappa*0.000001+2*alpha^2/(beta));
B_delta=(alpha^2/(beta));
delta_opt=min(A_delta/(2*B_delta),2);

%% Solver
counter=1;

%for  delta = [0.7:0.1:2.2,delta_opt]
    delta=1.5
L=alpha^2/((Kdr)*delta);

t=t_0;
Dirichletu(1:NxP2,2)=uyAnalytic(10,t);
f_10=@(x,y) f_1(x,y,t);
f_20=@(x,y) tau*f_2(x,y,t);
[p,Ap] = FEMParabolic2D(Coordinates,Elements,1/M+L,kappa*tau,Dirichletp,DirichletValue,Neumann,g,f_20,zeros(2*(2*Nx-1)*(2*Ny-1),1),CoordinatesP2,ElementsP2,1/M*p_0+L*p_0,1,0,Coefficients,GaussValues);
[u,Au] = FEMLinStrain2DvectorialP2(Coordinates,ElementsP2,2*mu,lambda,Dirichletu,Neumann,g,f_10,alpha*p,zeros(4*NN,1),1,0,CoefficientsP2,GaussValuesP2,CoordinatesP2);
u_old = u_0;
p_old = p_0;
errorp=norm(p-p_old,inf)/norm(p,inf);
erroru=norm(u-u_old,inf)/norm(u,inf);
iterations=1;
while errorp>10^(-6) || erroru>10^(-6)
    u_prev=u;
    p_prev=p;
    [p,~] = FEMParabolic2D(Coordinates,Elements,1/M+L,kappa*tau,Dirichletp,DirichletValue,Neumann,g,f_20,-alpha*(u_prev-u_old),CoordinatesP2,ElementsP2,1/M*p_old+L*p_prev,0,Ap,Coefficients,GaussValues);
    [u,~] = FEMLinStrain2DvectorialP2(Coordinates,ElementsP2,2*mu,lambda,Dirichletu,Neumann,g,f_10,alpha*p,zeros(4*NN,1),0,Au,CoefficientsP2,GaussValuesP2,CoordinatesP2);
    errorp=norm(p-p_prev,inf)/norm(p,inf)
    erroru=norm(u-u_prev,inf)/norm(u,inf)
    %error = sqrt(errorp^2+erroru^2)
    iterations=iterations+1;
end
t=t+tau;
while t<T+tau
    Dirichletu(1:NxP2,2)=uyAnalytic(10,t);
    u_old = u;
    p_old = p;
    f_1t = @(x,y) f_1(x,y,t);
    f_2t = @(x,y) tau*f_2(x,y,t)*scale;
    errorp = 1; erroru=1;
    while errorp>10^(-6) || erroru>10^(-6)
        u_prev=u;
        p_prev=p;
        [p,~] = FEMParabolic2D(Coordinates,Elements,1/M+L,kappa*tau,Dirichletp,DirichletValue,Neumann,g,f_20,-alpha*(u_prev-u_old),CoordinatesP2,ElementsP2,1/M*p_old+L*p_prev,0,Ap,Coefficients,GaussValues);
        [u,~] = FEMLinStrain2DvectorialP2(Coordinates,ElementsP2,2*mu,lambda,Dirichletu,Neumann,g,f_10,alpha*p,zeros(4*NN,1),0,Au,CoefficientsP2,GaussValuesP2,CoordinatesP2);
        errorp = norm(p-p_prev,inf)/norm(p,inf)
        erroru = norm(u-u_prev,inf)/norm(u,inf)
        %error = sqrt(errorp^2+erroru^2)
        iterations=iterations+1;
    end
    t=t+tau;
end
% Analysis(counter,2*index-1)=delta;
% Analysis(counter,2*index)=iterations;
% counter=counter+1;
% end
% 
% plot(Analysis(1:16,2*index-1),Analysis(1:16,2*index))
%  hold on
%  plot(Analysis(17,2*index-1),Analysis(17,2*index),'p','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',10)
%  end
%subplot(3,1,1)
%trisurf(Elements,X,Y,u(1:2:2*NN-1))
%subplot(3,1,2)
%trisurf(Elements,X,Y,u(2:2:2*NN))
%subplot(3,1,3)
%trisurf(Elements,X,Y,p)
usq=abs(u(1:2:2*NN-1).*u(2:2:2*NN)); 
subplot(1,2,1)
title('|u|')
trisurf(Elements,X,Y,usq,'EdgeColor','none','FaceColor','interp')
% subplot(3,2,1)
% trisurf(Elements,X,Y,u(2:2:2*NN))
subplot(1,2,2)
title('p')
trisurf(Elements,X,Y,p,'EdgeColor','none','FaceColor','interp')