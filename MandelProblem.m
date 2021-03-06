%% Mandel


%% Time
t_0=0;
tau=10^(0);
T=10^(1)*0.5;

%% Problem
a = 100;
b = 10; 
E =  5.94*10^9;% [Pa]
v = 0.2;% 
cf=3.03e-10;% [1/Pa]
alpha=1.0;% 
kappa=10^(-10);
B=0.833333; %Skempton coefficient 
vu= 0.44; %Undrained Poison's ratio
M=1.65*10^10;% Biot Modulus [Pa]
F=6*10^8;
K=lambda+2/3*mu;
Ku= K + alpha^2*M;
cf = M*kappa*(K+4/3*mu)/(Ku+4/3*mu);
lambda=(E*v)/((1-2*v)*(1+v));
mu= E/(2+2*v);
visc=0.01;% Viscocity

Beta=((1-2*v)*(1+v))/(2*E*(1-v));
kdr=1/Beta;
%mu=1;
%c=0.465; %[m2/s]
%k=100;
c=(2*kappa*(B^2)*mu*(1-v)*(1+vu)^2) /(9*visc*(1-vu)*(vu-v));
PresNormalizer=2*F*B*(1+vu)/(3*a);
%lambda = 1; mu=1; M=1; alpha=1; kappa=10^(0);
u_0=zeros(2*NN,1);
p_0=zeros(NN,1);
f_1=@(x,y,t) [0;0];
f_2=@(x,y,t) 0;

%% Analytical solution

[pAnalytic,uxAnalytic,uyAnalytic]=MandelAnalytic(F,B,vu,100,cf,mu,v,300);

%% MESH
x_min=0; x_max=100; y_min=0; y_max=10;
dx=2.5; dy=2.5;
[x,y]=meshgrid(x_min:dx:x_max,y_min:dy:y_max);
X=reshape(x,[],1); Y=reshape(y,[],1);
DT=delaunayTriangulation(X,Y);
Coordinates=DT.Points;
Elements=DT.ConnectivityList;
NN=length(X);
Neumann=[];
g=0;
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
%Boundary for flow
Dirichletp=find(Coordinates(:,1)==x_max);
DirichletValue = 0;
% Neumannp = zeros(3*(1/h+1),2);
% Neumannp(:,1)=[find(Coordinates(:,1)==x_min);find(Coordinates(:,2)==y_min);find(Coordinates(:,2)==y_max)];
% Neumannp(:,2)=0;
%Boundary for mechanics
Dirichletu = zeros(2*(100/dx+1)+10/dy+1,2);
Dirichletu(:,1)=[2*find(Coordinates(:,1)==x_min)-ones(length(find(Coordinates(:,1)==x_min)),1);2*find(Coordinates(:,2)==y_min);2*find(Coordinates(:,2)==y_max)];
Dirichletu(1:(10/dy+1),2)=0;
Dirichletu(10/dy+2:10/dy+1+(100/dx+1),2)=0;



%% Solver
t=t_0;
Dirichletu(10/dy+1+(100/dx+1)+1:10/dy+1+2*(100/dx+1),2)=uyAnalytic(10,t);
delta=2;
L=alpha^2/((mu+lambda)*delta);
f_10=@(x,y) f_1(x,y,t);
f_20=@(x,y) tau*f_2(x,y,t);
[p,Ap] = FEMParabolic2D(Coordinates,Elements,1/M+L,kappa*tau,Dirichletp,DirichletValue,Neumann,g,f_20,zeros(2*NN,1),1/M*p_0+L*p_0,1,0,Coefficients,GaussValues);
[u,Au] = FEMLinStrain2DvectorialP2(Coordinates,Elements,2*mu,lambda,Dirichletu,Neumann,g,f_10,alpha*p,zeros(2*NN,1),1,0,CoefficientsP2,GaussValuesP2);
u_old = u_0;
p_old = p_0;
errorp=norm(p-p_old,inf)/norm(p,inf);
erroru=norm(u-u_old,inf)/norm(u,inf);
iterations=1;
while errorp>10^(-12) || erroru>10^(-12)
    u_prev=u;
    p_prev=p;
    [p,~] = FEMParabolic2D(Coordinates,Elements,1/M+L,kappa*tau,Dirichletp,DirichletValue,Neumann,g,f_20,-alpha*(u_prev-u_old),1/M*p_old+L*p_prev,0,Ap,Coefficients,GaussValues);
    [u,~] = FEMLinStrain2DvectorialP2(Coordinates,Elements,2*mu,lambda,Dirichletu,Neumann,g,f_10,alpha*p,zeros(2*NN,1),0,Au,CoefficientsP2,GaussValuesP2);
    errorp=norm(p-p_prev,inf)/norm(p,inf)
    erroru=norm(u-u_prev,inf)/norm(u,inf)
    %error = sqrt(errorp^2+erroru^2)
    iterations=iterations+1;
end
t=t+tau;
while t<T+tau
    Dirichletu(10/dy+1+(100/dx+1)+1:10/dy+1+2*(100/dx+1),2)=uyAnalytic(10,t);
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
        errorp = norm(p-p_prev,inf)/norm(p,inf);
        erroru = norm(u-u_prev,inf)/norm(u,inf);
        %error = sqrt(errorp^2+erroru^2)
        iterations=iterations+1;
    end
    t=t+tau
end

subplot(4,1,1)
trisurf(Elements,X,Y,u(1:2:2*NN-1))
subplot(4,1,2)
trisurf(Elements,X,Y,u(2:2:2*NN))
subplot(4,1,3)
trisurf(Elements,X,Y,p/PresNormalizer)


% p=zeros(10000,1);
% for i=0:0.01:100
% p(i*100)=pAnalytic(i,T);
% end
p=[];
for i=0:0.01:100
    p=[p,pAnalytic(i,T)];
end
subplot(4,1,4)
plot(0:10^(-4):1,p/PresNormalizer);
    
