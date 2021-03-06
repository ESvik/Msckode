%%L - shape

%% Biot Analysis
%% MESH
x_min=0; x_max=1; y_min=0; y_max=1;
h=1/8;
[x1,y1]=meshgrid(x_min:h:x_max,y_min:h:y_max/2);
[x2,y2]=meshgrid(x_min:h:x_max/2,y_max/2:h:y_max);


X1=reshape(x1,[],1); Y1=reshape(y1,[],1);
X2=reshape(x2,[],1); Y2=reshape(y2,[],1);

[x,y]=meshgrid(x_min:h:x_max,y_min:h:y_max);
Xsynthetic=reshape(x,[],1);

% Xmid=find(X==0.5);
% Ymid=find(Y==0.5);
% Yxmid=find(Y(Xmid)==0.5);
% Ytop=find(Y(Xmid)==1);
% Xleft=find(X(Ymid)==1);
% top=Xmid(Ytop);
% mid=Xmid(Yxmid);
% left=Ymid(Xleft);
% Constraints=[mid,top;mid,left];
DT1=delaunayTriangulation(X1,Y1);
DT2=delaunayTriangulation(X2,Y2);


Coordinates1=DT1.Points;
Coordinates2=DT2.Points;
Elements1=DT1.ConnectivityList;
Elements3=DT2.ConnectivityList+length(Coordinates1);
Elements2=DT2.ConnectivityList+length(Coordinates1);

GoodPointsPre=find(Coordinates1(:,2)==y_max/2);
GoodPoints=zeros(length(x_min:h:x_max/2),1);
Goodline=x_min:h:x_max/2;
for i=1:length(x_min:h:x_max/2)
   GoodPoints(i)=GoodPointsPre(find(Coordinates1(GoodPointsPre,1)==Goodline(i)));
end

BadPoints=find(Coordinates2(:,2)==y_max/2);
Coordinates2(BadPoints,:)=[];
BadPoints=BadPoints+length(Coordinates1);
for i=1:(length(BadPoints)-1)
    for n=1:length(Elements2(:,1))
        if Elements2(n,1)==BadPoints(i)
            Elements2(n,1)=GoodPoints(i);
        end
        if Elements2(n,2)==BadPoints(i)
            Elements2(n,2)=GoodPoints(i);
        end
        if Elements2(n,3)==BadPoints(i)
            Elements2(n,3)=GoodPoints(i);
        end
        if Elements2(n,1)>BadPoints(i) && Elements2(n,1)<BadPoints(i+1)
            Elements2(n,1)=Elements2(n,1)-i;
        end
        if Elements2(n,2)>BadPoints(i) && Elements2(n,2)<BadPoints(i+1)
            Elements2(n,2)=Elements2(n,2)-i;
        end
        if Elements2(n,3)>BadPoints(i) && Elements2(n,3)<BadPoints(i+1)
            Elements2(n,3)=Elements2(n,3)-i;
        end
    end
end
for n=1:length(Elements2(:,1))
        if Elements2(n,1)==BadPoints(end)
            Elements2(n,1)=GoodPoints(end);
        end
        if Elements2(n,2)==BadPoints(end)
            Elements2(n,2)=GoodPoints(end);
        end
        if Elements2(n,3)==BadPoints(end)
            Elements2(n,3)=GoodPoints(end);
        end
        if Elements2(n,1)>BadPoints(end)
            Elements2(n,1)=Elements2(n,1)-length(BadPoints);
        end
        if Elements2(n,2)>BadPoints(end)
            Elements2(n,2)=Elements2(n,2)-length(BadPoints);
        end
        if Elements2(n,3)>BadPoints(end)
            Elements2(n,3)=Elements2(n,3)-length(BadPoints);
        end
end        
            

Elements=[Elements1;Elements2];
Coordinates=[Coordinates1;Coordinates2];
% Coordinates2(BadPoints,:)=[];
% Coordinates2=[Coordinates1;Coordinates2];
% Coordinates3=DT3.Points+length(Coordinates2);
% Coordinates=[Coordinates2;Coordinates3];
% Elements1=DT1.ConnectivityList;
% Elements2=DT2.ConnectivityList+length(Coordinates1);
% Elements3=DT3.ConnectivityList+length(Coordinates2);
% Elements=[Elements1;Elements2;Elements3];

Neumann=[];
g=0;
NN=length(Coordinates);
NNnew = length(Xsynthetic);
[ElementsP2,CoordinatesP2]=ElementsPlusEdgesLshape(Elements,Coordinates,NN,NNnew,x_max,h);
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
% Dirichletp = boundary(Coordinates(:,1),Coordinates(:,2),1);
% Dirichlet1=find(Coordinates(:,1)==0.5);
% Dirichletp=[Dirichletp;Dirichlet1(find(Coordinates(Dirichlet1,2)==0.5))];
% 
% 
% Dirichletupre=boundary(CoordinatesP2(:,1),CoordinatesP2(:,2),1);
% Dirichlet2=find(CoordinatesP2(:,1)==0.5);
% Dirichletupre=[Dirichletupre;Dirichlet2(find(CoordinatesP2(Dirichlet2,2)==0.5))];
% 
% Dirichletu = zeros(2*length(Dirichletupre),2);
% Dirichletu(1:2:end,1)=2*Dirichletupre-1;
% Dirichletu(2:2:end,1)=2*Dirichletupre;

Edgevector=0.5:h:1;
EdgevectorP2=0.5:h/2:1;
XCoordinates=find(Coordinates(:,1)==0.5);
YCoordinates=find(Coordinates(:,2)==0.5);
XCoordinatesP2=find(CoordinatesP2(:,1)==0.5);
YCoordinatesP2=find(CoordinatesP2(:,2)==0.5);
DirichletInnerEdgeUp=zeros(length(Edgevector),1);
DirichletInnerEdgeRight=zeros(length(Edgevector),1);
DirichletInnerEdgeUpP2=zeros(length(EdgevectorP2),1);
DirichletInnerEdgeRightP2=zeros(length(EdgevectorP2),1);
for i = 1:length(Edgevector)
    DirichletInnerEdgeRight(i)=XCoordinates(find(Coordinates(XCoordinates,2)==Edgevector(i)));
    DirichletInnerEdgeUp(i)=YCoordinates(find(Coordinates(YCoordinates,1)==Edgevector(i)));
end
for i =1:length(EdgevectorP2)
DirichletInnerEdgeRightP2(i)=XCoordinatesP2(find(CoordinatesP2(XCoordinatesP2,2)==EdgevectorP2(i)));
DirichletInnerEdgeUpP2(i)=YCoordinatesP2(find(CoordinatesP2(YCoordinatesP2,1)==EdgevectorP2(i)));
end
DirichletTop=find(Coordinates(:,2)==y_max);
DirichletLeft=find(Coordinates(:,1)==x_min);
DirichletDown=find(Coordinates(:,2)==y_min);
DirichletRight=find(Coordinates(:,1)==x_max);
%DirichletTopP2=find(CoordinatesP2(:,2)==y_max);
DirichletLeftP2=find(CoordinatesP2(:,1)==x_min);
DirichletDownP2=find(CoordinatesP2(:,2)==y_min);
DirichletRightP2=find(CoordinatesP2(:,1)==x_max);

Dirichletupre=[DirichletDownP2;DirichletRightP2;DirichletInnerEdgeUpP2;DirichletInnerEdgeRightP2;DirichletLeftP2];
Dirichletupre=[2*Dirichletupre;2*Dirichletupre-1];
Dirichletu=zeros(length(Dirichletupre),2);
Dirichletu(:,1)=Dirichletupre;
Dirichletp=[DirichletDown;DirichletRight;DirichletInnerEdgeUp;DirichletInnerEdgeRight;DirichletTop;DirichletLeft];

DirichletValue=0;

%% Time
t_0=0;
tau=10^(-1);
T=tau;

%% Problem
kappavector = [10^(-15),10^(-14),10^(-13),10^(-12),10^(-11),10^(-10)];
Analysis=zeros(40,12);
%for index=1:6
kappa = kappavector(6);
%pressurescale=1/kappa*10^(-4);
pressurescale = 10^11;
uexact = @(x,y,t) [t.*x.*y.*(x-1).*(y-1),t.*x.*y.*(x-1).*(y-1)];
pexact = @(x,y,t) pressurescale*t.*x.*y.*(x-1).*(y-1);
lambda = 27.778*10^(9); mu=41.667*10^(9); M=100*10^9; alpha=1; %kappa=10^(-10);
%lambda = 1; mu=1; M=1; alpha=1; kappa=10^(0);
u_0=zeros(2*length(CoordinatesP2(:,1)),1);
u0=uexact(CoordinatesP2(:,1),CoordinatesP2(:,2),t_0);
u_0(1:2:length(u_0))=u0(:,1);
u_0(2:2:length(u_0))=u0(:,2);
p_0=pexact(Coordinates(:,1),Coordinates(:,2),t_0);
%f_1=@(x,y,t) [0;0];
%f_2=@(x,y,t) 0;
f_1=@(x,y,t) [(-2*mu-lambda)*2*t*y.*(y-1)+(-mu-lambda)*(2*x-1).*(2*y-1).*t-mu*2*t*x.*(x-1)+(alpha*t*y.*(y-1).*(2*x-1))*pressurescale; -mu*2*t*y.*(y-1)+(-mu-lambda)*t*(2*x-1).*(2*y-1)+(-2*mu-lambda)*2*t*x.*(x-1)+(alpha*t*x.*(x-1).*(2*y-1))*pressurescale];
f_2=@(x,y,t) (1/M*x.*y.*(x-1).*(y-1)-kappa*t*2*(x.*(x-1)+y.*(y-1)))*pressurescale+alpha*(y.*(y-1).*(2*x-1)+x.*(x-1).*(2*y-1));

%% Mathematical optima

Kdr=1.4*mu+lambda;

beta = Kdr;
A_delta=(2/M+2*tau*kappa+2*alpha^2/(beta));
B_delta=(alpha^2/(beta));
delta_opt=min(A_delta/(2*B_delta),2);


%% Solver
counter=1;
%for  delta = [0.7:0.05:2.6,delta_opt]
    delta=1.5
    %L=alpha^2/((mu+lambda)*delta);
    L=alpha^2/(Kdr*delta);
    t=t_0+tau; 
    f_10=@(x,y) f_1(x,y,t);
    f_20=@(x,y) tau*f_2(x,y,t);
    [p,Ap] = FEMParabolic2D(Coordinates,Elements,1/M+L,kappa*tau,Dirichletp,DirichletValue,Neumann,g,f_20,zeros(2*length(CoordinatesP2(:,1)),1),CoordinatesP2,ElementsP2,1/M*p_0+L*p_0,1,0,Coefficients,GaussValues);
    [u,Au] = FEMLinStrain2DvectorialP2(Coordinates,ElementsP2,2*mu,lambda,Dirichletu,Neumann,g,f_10,alpha*p,zeros(4*NN,1),1,0,CoefficientsP2,GaussValuesP2,CoordinatesP2);
    u_old = u_0;
    p_old = p_0;
    errorp=norm(p-p_old,inf)/norm(p,inf);
    erroru=norm(u-u_old,inf)/norm(u,inf);
    iterations=1;
    while errorp>10^(-12) || erroru>10^(-12)
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
    u_old = u;
    p_old = p;
    f_1t = @(x,y) f_1(x,y,t);
    f_2t = @(x,y) tau*f_2(x,y,t);
    errorp = 1; erroru=1;
    while errorp>10^(-12) || erroru>10^(-12)
        u_prev=u;
        p_prev=p;
        [p,~] = FEMParabolic2D(Coordinates,Elements,1/M+L,kappa*tau,Dirichletp,DirichletValue,Neumann,g,f_20,-alpha*(u_prev-u_old),CoordinatesP2,ElementsP2,1/M*p_old+L*p_prev,0,Ap,Coefficients,GaussValues);
        [u,~] = FEMLinStrain2DvectorialP2(Coordinates,ElementsP2,2*mu,lambda,Dirichletu,Neumann,g,f_10,alpha*p,zeros(4*NN,1),0,Au,CoefficientsP2,GaussValuesP2,CoordinatesP2);
        errorp = norm(p-p_prev,inf)/norm(p,inf);
        erroru = norm(u-u_prev,inf)/norm(u,inf);
        %error = sqrt(errorp^2+erroru^2)
        iterations=iterations+1;
    end
    t=t+tau;
end
% Analysis(counter,2*index-1)=delta;
% Analysis(counter,2*index)=iterations;
% counter=counter+1;
% end
% plot(Analysis(1:39,2*index-1),Analysis(1:39,2*index))
% hold on
% plot(Analysis(40,2*index-1),Analysis(40,2*index),'p','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',10)
% end

%subplot(3,1,1)
%trisurf(Elements,Coordinates(:,1),Coordinates(:,2),u(1:2:2*NN-1))
%subplot(3,1,2)
%trisurf(Elements,Coordinates(:,1),Coordinates(:,2),u(2:2:2*NN))
%subplot(3,1,3)
%trisurf(Elements,Coordinates(:,1),Coordinates(:,2),p)
usq=abs(u(1:2:2*NN-1).*u(2:2:2*NN)); 
subplot(1,2,1)
title('|u|')
trisurf(Elements,Coordinates(:,1),Coordinates(:,2),usq,'EdgeColor','none','FaceColor','interp')
% subplot(3,2,1)
% trisurf(Elements,X,Y,u(2:2:2*NN))
subplot(1,2,2)
title('p')
trisurf(Elements,Coordinates(:,1),Coordinates(:,2),p,'EdgeColor','none','FaceColor','interp')

