%% 2D vectorial test
%% Mesh
x_min=0; x_max=1; y_min=0; y_max=1;
[x,y]=meshgrid(x_min:0.05:x_max,y_min:0.05:y_max);
X=reshape(x,[],1); Y=reshape(y,[],1);
DT=delaunayTriangulation(X,Y);
Coordinates=DT.Points;
Elements=DT.ConnectivityList;
Neumann=[];
Dirichlet=[2*find(Coordinates(:,1)==x_min);2*find(Coordinates(:,1)==x_min)-ones(length(find(Coordinates(:,1)==x_min)),1);2*find(Coordinates(:,1)==x_max);2*find(Coordinates(:,1)==x_max)-ones(length(find(Coordinates(:,1)==x_max)),1);2*find(Coordinates(:,2)==y_min);2*find(Coordinates(:,2)==y_min)-ones(length(find(Coordinates(:,2)==y_min)),1);2*find(Coordinates(:,2)==y_max);2*find(Coordinates(:,2)==y_max)-ones(length(find(Coordinates(:,2)==y_max)),1)];
DirichletValue=0;
NN=length(Coordinates(:,1));
uexact=@(x,y) x.*(x-1).*y.*(y-1);
f=@(x,y) [x.*(1-x).*y.*(1-y)+2*y.*(1-y)-(2*x-1).*(2*y-1);x.*(1-x).*y.*(1-y)+2*x.*(1-x)-(2*x-1).*(2*y-1)];
%f=@(x,y) [2*y.*(1-y)-(1-y).*(1-x)-x.*y+x.*(1-y)+y.*(1-x);2*y.*(1-y)-(1-y).*(1-x)-x.*y+x.*(1-y)+y.*(1-x);];
%f=@(x,y) [x.*(x-1).*y.*(y-1)-2*y.^2+4*y+2*x-4*x.*y-1;x.*(x-1).*(y-1).*y-2*x.^2+4*x+2*y-4*x.*y-1];
%f = @(x,y) [x.*(x-1).*y.*(y-1)-2*x.*(x-1)-2*y.*(y-1);x.*(x-1).*y.*(y-1)-2*x.*(x-1)-2*y.*(y-1)];
u=FEMParabolic2Dvectorial(Coordinates,Elements,1,1,Dirichlet,DirichletValue,Neumann,0,f,0,zeros(2*NN,1));
%Wolframalpha av f-u:(1 - 4 y + 2 y^2 + x (-2 + 4 y), 1 + 2 x^2 + 4 x (-1 + y) - 2 y) 
error1=norm(u(1:2:2*NN-1)-uexact(X,Y),inf);
error2=norm(u(2:2:2*NN)-uexact(X,Y),inf);

subplot(2,1,1)
trisurf(Elements,X,Y,u(1:2:2*NN-1))
hold on
subplot(2,1,2)
trisurf(Elements,X,Y,u(1:2:2*NN-1))

