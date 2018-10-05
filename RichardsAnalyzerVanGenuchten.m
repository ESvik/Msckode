%Richards Eq in 2D Nonlinear Diffusion

%% Mesh
x_min=0; x_max=1; y_min=0; y_max=1;

h=1/4;
[x,y]=meshgrid(x_min:h:x_max,y_min:h:y_max);
X=reshape(x,[],1); Y=reshape(y,[],1);
DT=delaunayTriangulation(X,Y);
Coordinates=DT.Points;
Elements=DT.ConnectivityList;
Neumann=[];
Dirichlet=zeros(3*(1/h+1),2);
Dirichlet(:,1)=[find(Coordinates(:,1)==x_min);find(Coordinates(:,1)==x_max);find(Coordinates(:,2)==y_min)];
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
    GaussValues(i,1) = Coefficients(i,1) + Coefficients(i,2)*1/6 + Coefficients(i,3)*1/6;
    GaussValues(i,2) = Coefficients(i,1) + Coefficients(i,2)*2/3 + Coefficients(i,3)*1/6;
    GaussValues(i,3) = Coefficients(i,1) + Coefficients(i,2)*1/6 + Coefficients(i,3)*2/3;
end

%% Time
t_0=0;
tau=0.1;
T=t_0+tau;

%% Problem

%Coefficients


kappaabs =3*10^(-2); mu=1; a=0.1844; n=3;
k=@(s) kappaabs/mu*sqrt(s).*(1-(1-s.^(n/(n-1))).^((n-1)/n)).^2;
b=@(u) genuchtensaturation(u,a,n);
kappa=@(u) k(b(u));
%b = @(u) step(u,steps,L_b,b_m); 
bprime = @(u) genuchtensaturationprime(u,a,n);
bdoubleprime = @(u) genuchtensaturationdoubleprime(u,a,n);
%kappaprime = @(u) (kappaabs/mu*0.5*b(u).^(-1/2).*(1-(1-b(u).^(n/n-1)).^((n-1)/n)).^2+2*kappaabs/mu*sqrt(b(u)).*(1-(1-b(u).^(n/(n-1))).^((n-1)/n)).*(1-b(u).^(n/(n-1))).^((n-1)/n-1).*b(u).^(n/(n-1)-1)).*bprime(u);
f1= @(x,y,t) 0;
f2= @(x,y,t) f1(x,y,t)*tau;
g=0;

bprimemax=max(bprime(-10:0.1:0));
bm=min(bprime(-10:0.01:0.1));

L_b=bprimemax;
%kappaprime =@(s) ((kappaabs*(1 - (1 - s.^(n/(n - 1))).^((n - 1)/n)).^2)./(2*mu.*sqrt(s)) + (2*kappaabs*s.^(n/(n - 1) - 1/2).*(1 - (1 - s.^(n/(n - 1))).^((n - 1)/n)).*(1 - s.^(n/(n - 1))).^((n - 1)/n - 1))./mu);
epsilon=0.01;
kappaprime=@(u) (kappa(u)-kappa(u-epsilon))/epsilon;

%% Solver
KIRCHHOFF = 2;
Errors=zeros(5,1);
for SCHEME=1:4

% L=L_b;
% 
% t=t_0+tau;
% u_0 = -7.78*ones(length(X),1);
% f=@(x,y) f2(x,y,t);

% [u,Au]=FEMParabolic2DP1Richards(Coordinates,Elements,L,0,Dirichlet,Neumann,g,f,b(u_0)-b(u_0),1,0,Coefficients,GaussValues,bprime,u_0,SCHEME,0,KIRCHHOFF,kappa,kappaprime,tau);
% u_old = u_0;
% error=norm(u-u_old,inf)+norm(u-u_old,inf)/norm(u,inf);
% %Errors(hhh,1)=log(error);
% errorcounter=2;
% iterations=1;
% while error > 10^(-8)
%         u_prev=u;
%         LAMBDA=norm(u_prev-u_old,inf)/tau;
%         ModifiedL=LAMBDA*tau*norm(bdoubleprime(u_prev),inf)*0.1;
%         [u,Au]=FEMParabolic2DP1Richards(Coordinates,Elements,L,0,Dirichlet,Neumann,g,f,b(u_old)-b(u_prev),1,0,Coefficients,GaussValues,bprime,u_prev,SCHEME,ModifiedL,KIRCHHOFF,kappa,kappaprime,tau);
%         error=norm(u-u_prev,inf)+norm(u-u_prev,inf)/norm(u,inf)
%        % Errors(hhh,errorcounter)=log(error);
%         errorcounter=errorcounter+1;
%         iterations=iterations+1;
% end
% 
% uPre=u;
  
L=0.5*L_b;

t=t_0+tau;
u_0 = -7.78*ones(length(X),1);
f=@(x,y) f2(x,y,t);

[u,Au]=FEMParabolic2DP1Richards(Coordinates,Elements,L,0,Dirichlet,Neumann,g,f,b(u_0)-b(u_0),1,0,Coefficients,GaussValues,bprime,u_0,SCHEME,0,KIRCHHOFF,kappa,kappaprime,tau);
u_old = u_0;
error=norm(u-u_old,inf)+norm(u-u_old,inf)/norm(u,inf);
Errors(SCHEME,1)=log(error);
errorcounter=2;
iterations=1;
while error > 10^(-8)
        u_prev=u;
        LAMBDA=norm(u_prev-u_old,inf)/tau;
        ModifiedL=LAMBDA*tau*norm(bdoubleprime(u_prev),inf)*0.1;
        [u,Au]=FEMParabolic2DP1Richards(Coordinates,Elements,L,0,Dirichlet,Neumann,g,f,b(u_old)-b(u_prev),1,0,Coefficients,GaussValues,bprime,u_prev,SCHEME,ModifiedL,KIRCHHOFF,kappa,kappaprime,tau);
        error=norm(u-u_prev,inf)+norm(u-u_prev,inf)/norm(u,inf)
        Errors(SCHEME,errorcounter)=log(error);
        errorcounter=errorcounter+1;
        iterations=iterations+1;
end

end

 %Picard
% t=t_0+tau;
% u_0 = -7.78*ones(length(X),1);
% f=@(x,y) f2(x,y,t);
% L=bprimemax;
% [u,Au]=FEMParabolic2DP1Richards(Coordinates,Elements,L,0,Dirichlet,Neumann,g,f,b(u_0)-b(u_0),1,0,Coefficients,GaussValues,bprime,u_0,2,0,KIRCHHOFF,kappa,kappaprime,tau);
% u_old = u_0;
% error=norm(u-u_old,inf)+norm(u-u_old,inf)/norm(u,inf);
% errorcounter=2;
% iterations=1;
% while error > 10^(-8)
%         u_prev=u;
%         [u,Au]=FEMParabolic2DP1Richards(Coordinates,Elements,L,0,Dirichlet,Neumann,g,f,b(u_old)-b(u_prev),0,Au,Coefficients,GaussValues,bprime,u_prev,2,0,KIRCHHOFF,kappa,kappaprime,tau);
%         error=norm(u-u_prev,inf)+norm(u-u_prev,inf)/norm(u,inf);
%         errorcounter=errorcounter+1;
%         iterations=iterations+1;
% end
% 
% uPre=u;
% 
t=t_0+tau;
u_0 = -7.78*ones(length(X),1);
f=@(x,y) f2(x,y,t);
L=bprimemax;
[u,Au]=FEMParabolic2DP1Richards(Coordinates,Elements,L,0,Dirichlet,Neumann,g,f,b(u_0)-b(u_0),1,0,Coefficients,GaussValues,bprime,u_0,2,0,KIRCHHOFF,kappa,kappaprime,tau);
u_old = u_0;
error=norm(u-u_old,inf)+norm(u-u_old,inf)/norm(u,inf);
Errors(5,1)=log(error);
errorcounter=2;
iterations=1;
while error > 10^(-8)
        u_prev=u;
        [u,Au]=FEMParabolic2DP1Richards(Coordinates,Elements,L,0,Dirichlet,Neumann,g,f,b(u_old)-b(u_prev),0,Au,Coefficients,GaussValues,bprime,u_prev,2,0,KIRCHHOFF,kappa,kappaprime,tau);
        error=norm(u-u_prev,inf)+norm(u-u_prev,inf)/norm(u,inf);
        Errors(5,errorcounter)=log(error);
        errorcounter=errorcounter+1;
        iterations=iterations+1;
end

% 
% order = zeros(5,1);
% for schemes=1:5
% for i =3:length(nonzeros(Errors(schemes,:)))
%     order(schemes,i-2)=log(Errors(schemes,i)/Errors(schemes,i-1))/log(Errors(schemes,i-1)/Errors(schemes,i-2));
% end
% end
% p1=plot(1+2:2+length(nonzeros(order(1,:))),nonzeros(order(1,:)),'b-o');
% hold on
% p2=plot(1+2:2+length(nonzeros(order(2,:))),nonzeros(order(2,:)),'r-.');
% p3=plot(1+2:2+length(nonzeros(order(3,:))),nonzeros(order(3,:)),'g-*');
% p4=plot(1+2:2+length(nonzeros(order(4,:))),nonzeros(order(4,:)),'c-p');
% p5=plot(1+2:2+length(nonzeros(order(5,:))),nonzeros(order(5,:)),'k--');
% 
% legend([p1,p2,p3,p4,p5],'L=L_s','Modified L','Newton','Locally optimized L','Modified Picard')
% ylabel('Order')
% xlabel('Iteration i')
% set(gca,'fontsize',25)

p1=plot(1:length(nonzeros(Errors(1,:))),nonzeros(Errors(1,:)),'b-.')
hold on
p2=plot(1:length(nonzeros(Errors(2,:))),nonzeros(Errors(2,:)),'r-*')
p3=plot(1:length(nonzeros(Errors(3,:))),nonzeros(Errors(3,:)),'g-o')
p4=plot(1:length(nonzeros(Errors(4,:))),nonzeros(Errors(4,:)),'c-p')
p5=plot(1:length(nonzeros(Errors(5,:))),nonzeros(Errors(5,:)),'k--')
legend([p1,p2,p3,p4,p5],'L=0.5L_s','Modified L','Newton','Locally optimized L','Modified Picard')
ylabel('log_{10}(||u^{i}-u^{i-1}||)')
xlabel('Iteration i')
set(gca,'fontsize',25)




% while t < T+tau
%     f= @(x,y) f2(x,y,t);
%     u_old=u;
%     error = 1;
%     while error > 10^(-8)
%         u_prev=u;
%         [u,Au]=FEMParabolic2DP1Richards(Coordinates,Elements,kappa,tau,Dirichlet,DirichletValue,Neumann,g,f,b(u_old)-b(u_prev),1,Au,Coefficients,GaussValues,bprime,u_prev);
%         error=norm(u-u_prev,inf)+norm(u-u_prev,inf)/norm(u,inf);
%         iterations=iterations+1;
%     end
%     t=t+tau;
% end
% analysis(counter,2*index-1)=1.1; analysis(counter,2*index)=iterations;
% counter=counter+1;
% gammaopt=max(0,1+(L_b*b_m-sqrt(L_b^2*b_m^2+4*L_b*b_m^2*kappa*tau+4*L_b*b_m*kappa^2*tau^2))/(4*b_m*kappa*tau));
% 
% for gamma=[0:0.05:0.5,gammaopt]
%     kappa
%     gamma
% L=L_b/(2*(1-gamma));
% t=t_0+tau;
% u_0 = uexact(X,Y,t_0);
% f=@(x,y) f2(x,y,t);
% [u,Au]=FEMParabolic2DP1Richards(Coordinates,Elements,L,kappa*tau,Dirichlet,DirichletValue,Neumann,g,f,b(u_0)-b(u_0)+L*u_prev,1,0,Coefficients,GaussValues);
% u_old = u_0;
% error=norm(u-u_old,inf)+norm(u-u_old,inf)/norm(u,inf);
% iterations=1;
% while error > 10^(-8)
%         u_prev=u;
%         [u,Au]=FEMParabolic2DP1Richards(Coordinates,Elements,L,kappa*tau,Dirichlet,DirichletValue,Neumann,g,f,b(u_old)-b(u_prev)+L*u_prev,0,Au,Coefficients,GaussValues);
%         error=norm(u-u_prev,inf)+norm(u-u_prev,inf)/norm(u,inf);
%         iterations=iterations+1;
% end
% t=t+tau;
% 
% while t < T+tau
%     f= @(x,y) f2(x,y,t);
%     u_old=u;
%     error = 1;
%     while error > 10^(-8)
%         u_prev=u;
%         [u,Au]=FEMParabolic2DP1Richards(Coordinates,Elements,l,kappa*tau,Dirichlet,DirichletValue,Neumann,g,f,b(u_old)-b(u_prev)+L*u_prev,0,Au,Coefficients,GaussValues);
%         error=norm(u-u_prev,inf)+norm(u-u_prev,inf)/norm(u,inf);
%         iterations=iterations+1;
%     end
%     t=t+tau;
% end
% analysis(counter,2*index-1)=L/L_b; analysis(counter,2*index)=iterations;
% counter=counter+1;
% end
% 
% end
% plot(Analysis(1,1:18)/L_b,Analysis(2,1:18)','k:')
% hold on
% plot(Analysis(1,19)/L_b,Analysis(2,19),'p','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',10)
% plot(Analysis(1,1:15)/L_b,Analysis(2,1:15)',':')
% hold on
%plot(Analysis(1,16)/L_b,Analysis(2,16),'p','MarkerEdgeColor','k','MarkerSize',10)
% legend('b_m=0.01','b_m=0.1','b_m=0.5','b_m=0.7','b_m=0.9')
%trisurf(Elements,X,Y,u)