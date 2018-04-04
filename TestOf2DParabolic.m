%2D Richards Eq

x_min=0; x_max=1; y_min=0; y_max=1;
[x,y]=meshgrid(x_min:0.01:x_max,y_min:0.01:y_max);
X=reshape(x,[],1); Y=reshape(y,[],1);
DT=delaunayTriangulation(X,Y);
Coordinates=DT.Points;
Elements=DT.ConnectivityList;
Neumann=[];
Dirichlet=[find(Coordinates(:,1)==x_min);find(Coordinates(:,1)==x_max);find(Coordinates(:,2)==y_min);find(Coordinates(:,2)==y_max)];
DirichletValue=0;
C=1; B=1;
f=@(x,y) ((1-x).*(1-y).*x.*y+2*((1-y).*y+(1-x).*x));
g=0;
fVector=zeros(length(Coordinates(:,1)),1);
uexact=@(x,y) x.*(x-1).*y.*(y-1);


u = FEMParabolic2D(Coordinates,Elements,C,B,Dirichlet,DirichletValue,Neumann,g,f,fVector);

error=norm(u-uexact(X,Y),2)

trisurf(Elements,X,Y,u)