%% Mesh
x_min=0; x_max=1; y_min=0; y_max=1;
[x,y]=meshgrid(x_min:0.05:x_max,y_min:0.05:y_max);
X=reshape(x,[],1); Y=reshape(y,[],1);
DT=delaunayTriangulation(X,Y);
Coordinates=DT.Points;
Elements=DT.ConnectivityList;
Neumann=[];
g=0;
Dirichlet=[find(Coordinates(:,1)==x_min);find(Coordinates(:,1)==x_max);find(Coordinates(:,2)==y_min);find(Coordinates(:,2)==y_max)];
DirichletValue=0;

%% Time
t_0=0.9;
tau=0.1;
T=1;

%% Problem
uexact = @(x,y,t) t.*x.*y.*(x-1).*(y-1);
vexact = @(x,y,t) -t.*x.*y.*(x-1).*(y-1);
lambda = 1; M=1; alpha=1; D_1=1; D_2=1;
u_0=uexact(X,Y,t_0);
v_0=vexact(X,Y,t_0);
f_1 = @(x,y,t) lambda*x.*t.*(x-1).*y.*(y-1)-2*D_1*t.*(y.*(y-1)+x.*(x-1))+alpha*x.*t.*(x-1).*y.*(y-1);
f_2 = @(x,y,t) -1/M*x.*(x-1).*y.*(y-1)+2*D_2*t.*(y.*(y-1)+x.*(x-1))+alpha*x.*(x-1).*y.*(y-1);


%%Solver
t=t_0+tau;
u=u_0;
v=v_0;
L=1;
while t<T+tau
    u_old = u;
    v_old = v;
    f_1 = @(x,y) f_1(x,y,t);
    f_2 = @(x,y) tau*f_2(x,y,t);
    error = 1;
    while error>0.00001
        u_prev=u;
        v_prev=v;
        v=FEMParabolic2D(Coordinates,Elements,1/M+L,D_2*tau,Dirichlet,DirichletValue,Neumann,g,f_2,1/M*v_old+alpha*u_old+L*v_prev-alpha*u_prev);
        u=FEMParabolic2D(Coordinates,Elements,lambda,D_1,Dirichlet,DirichletValue,Neumann,g,f_1,alpha*v);
        errorv=norm(v-v_prev,2); erroru=norm(u-u_prev,2);
        error = sqrt(errorv^2+erroru^2);
    end
    t=t+tau;
end

subplot(2,1,1)
trisurf(Elements,X,Y,u)
subplot(2,1,2)
trisurf(Elements,X,Y,v)
