%Richards Eq in 2D

%% Mesh
x_min=0; x_max=1; y_min=0; y_max=1;
[x,y]=meshgrid(x_min:0.05:x_max,y_min:0.05:y_max);
X=reshape(x,[],1); Y=reshape(y,[],1);
DT=delaunayTriangulation(X,Y);
Coordinates=DT.Points;
Elements=DT.ConnectivityList;
Neumann=[];
Dirichlet=[find(Coordinates(:,1)==x_min);find(Coordinates(:,1)==x_max);find(Coordinates(:,2)==y_min);find(Coordinates(:,2)==y_max)];
DirichletValue=0;
NN=length(Coordinates(:,1));

%% Time
t_0=7.9;
tau=0.1;
T=8;

%% Problem
uexact = @(x,y,t) x.*y.*(1-x).*(1-y).*t;
%Coefficients
k=0.1; L_b=10; steps=1; b_m=1;
pol=@(u) richardnonlin(u);
b= @(u) b_m.*u + (L_b-b_m)*(20/27)*pol(u);
%b = @(u) step(u,steps,L_b,b_m); 
f1= @(x,y,t) (b_m+(L_b-b_m)*(20/27).*(-5.4*uexact(x,y,t).^2+5.4*uexact(x,y,t))).*uexact(x,y,1)+2*k*t.*((1-y).*y+(1-x).*x);
f2= @(x,y,t) f1(x,y,t)*tau;
u_0 = @(x,y) 0;
g=0;
LL = @(g) 1./(2*(1-g));


%% Solver

Theoretical_Minima = theoreticalminimarichards(L_b,b_m,1,tau,k);
%if Theoretical_Minima > 0
L_Optb = LL(Theoretical_Minima);
%else
%    L_Optb =0.5;
%end
Analysis=zeros(2,8);
counter=1;
L=L_b*1
%    L/L_b
counter=1;
t=t_0+tau;
u = uexact(X,Y,t_0);
iterations=0;
    errorcomparison=[];
    errorfinished = zeros(6,1);
while t < T+tau
    f= @(x,y) f2(x,y,t);
    u_old=u;
    error = 1;
    while error > 0.00001
        u_prev=u;
        u=FEMParabolic2D(Coordinates,Elements,L,k*tau,Dirichlet,DirichletValue,Neumann,g,f,b(u_old)-b(u_prev)+L*u_prev);
        error=norm(u-u_prev,2);
        
%         errorn = 0;
%         for k = 1:length(Elements(:,1))
%             nodes = Elements(k,:);
%             vertices = Coordinates(nodes,:);
%             [ARefTri,bRefTri]=RefTriangleMapInv(vertices(:,1),vertices(:,2));
%             ARefTriInv = ARefTri^(-1);
%             VectCoeff=[1,vertices(1,1),vertices(1,2);1,vertices(2,1),vertices(2,2 );1,vertices(3,1),vertices(3,2)]\[abs(u(nodes(1))-u_prev(nodes(1))); abs(u(nodes(2))-u_prev(nodes(2)));abs(u(nodes(3))-u_prev(nodes(3)))];
%             GaussNodes= [ARefTri*[1/6;1/6]+bRefTri,ARefTri*[2/3;1/6]+bRefTri,ARefTri*[1/6;2/3]+bRefTri];
%             fVectGauss = [VectCoeff'*[1;GaussNodes(:,1)],VectCoeff'*[1;GaussNodes(:,2)],VectCoeff'*[1;GaussNodes(:,3)]];
%             errorn = errorn + det(ARefTri)*sum(1/6*fVectGauss);
%         end
%         error=errorn
    iterations = iterations+1;
    end
    t=t+tau;
end
uprecomp=u;
u = uexact(X,Y,t_0);




%end

% plot(Analysis(1,1:18)/L_b,Analysis(2,1:18)','k:')
% hold on
% plot(Analysis(1,19)/L_b,Analysis(2,19),'p','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',10)

% plot(Analysis(1,1:15)/L_b,Analysis(2,1:15)',':')
% hold on
%plot(Analysis(1,16)/L_b,Analysis(2,16),'p','MarkerEdgeColor','k','MarkerSize',10)


% legend('b_m=0.01','b_m=0.1','b_m=0.5','b_m=0.7','b_m=0.9')
%trisurf(Elements,X,Y,u)
 trisurf(Elements,X,Y,uexact(X,Y,T))