%% Linear Strain equation solver  C(e(u):e(v))+B(nabla*u,nabla*v)=(f,v) +(h,nabla * v)   in 2D domain with 2D Range using P2 elements.
function u = FEMLinStrain2DvectorialP2(Coordinates,Elements,C,B,Dirichlet,DirichletValue,Neumann,g,f,h,fVector)
    NN=length(Coordinates(:,1));
    A=sparse(2*NN,2*NN);
    b=zeros(2*NN,1);
    number_of_elements=length(Elements(:,1));
    
    %Basis and its evaluation points
    Coefficients= zeros(3,6);
    MatOfCoord=[1,0,0,0,0,0;1,1,0,1,0,0;1,0,1,0,1,0;1,1/2,0,1/4,0,0;1,0,1/2,0,1/4,0;1,1/2,1/2,1/4,1/4,1/4];
    Coefficients(1,:)=MatOfCoord\[1;0;0;1/2;1/2;0];
    Coefficients(2,:)=MatOfCoord\[0;1;0;1/2;0;1/2];
    Coefficients(3,:)=MatOfCoord\[0;0;1;0;1/2;1/2];
    
    Coefficients=[Coefficients(1,:);Coefficients(1,:);Coefficients(2,:);Coefficients(2,:);Coefficients(3,:);Coefficients(3,:)];
    GaussValues = zeros(6,3);
    for i = 1:6
        GaussValues(i,1) = (Coefficients(i,1) + Coefficients(i,2)*1/6 + Coefficients(i,3)*1/6+Coefficients(i,4)*1/6^2+Coefficients(i,5)*1/6^2+Coefficients(i,6)*1/6^2);
        GaussValues(i,2) = (Coefficients(i,1) + Coefficients(i,2)*2/3 + Coefficients(i,3)*1/6+Coefficients(i,4)*(2/3)^2+Coefficients(i,5)*1/6^2+Coefficients(i,6)*1/6*2/3);
        GaussValues(i,3) = (Coefficients(i,1) + Coefficients(i,2)*1/6 + Coefficients(i,3)*2/3+Coefficients(i,4)*1/6^2+Coefficients(i,5)*(2/3)^2+Coefficients(i,6)*1/6*2/3);
    end

    %Assembly
    for k = 1:number_of_elements
    nodes = Elements(k,:);
    vertices = Coordinates(nodes,:);
    hVect=[1,vertices(1,1),vertices(1,2);1,vertices(2,1),vertices(2,2);1,vertices(3,1),vertices(3,2)]\h(nodes);
    nodes=[2*nodes(1)-1,2*nodes(1),2*nodes(2)-1,2*nodes(2),2*nodes(3)-1,2*nodes(3)];
    
    
    %Transformation to refTriangle
    [ARefTri,bRefTri]=RefTriangleMapInv(vertices(:,1),vertices(:,2));
    [ARefTriinv,~]=RefTriangleMap(vertices(:,1),vertices(:,2));
    %ARefTriinv=ARefTri^(-1);
    
    %A
    ALdiv=zeros(6,6);
    for i =1:2:5
        for j = 1:2:5
             ALdiv(i,j)=1/6*B*det(ARefTri)*sum([(ARefTriinv(1,1)*(Coefficients(i,2)+2*Coefficients(i,4)*1/6+Coefficients(i,6)*1/6)+ARefTriinv(2,1)*(Coefficients(i,3)+2*Coefficients(i,5)*1/6+Coefficients(i,6)*1/6))*(ARefTriinv(1,1)*(Coefficients(j,2)+2*Coefficients(j,4)*1/6+Coefficients(j,6)*1/6)+ARefTriinv(2,1)*(Coefficients(j,3)+2*Coefficients(j,5)*1/6+Coefficients(j,6)*1/6)),(ARefTriinv(1,1)*(Coefficients(i,2)+2*Coefficients(i,4)*1/6+Coefficients(i,6)*2/3)+ARefTriinv(2,1)*(Coefficients(i,3)+2*Coefficients(i,5)*2/3+Coefficients(i,6)*1/6))*(ARefTriinv(1,1)*(Coefficients(j,2)+2*Coefficients(j,4)*1/6+Coefficients(j,6)*2/3)+ARefTriinv(2,1)*(Coefficients(j,3)+2*Coefficients(j,5)*2/3+Coefficients(j,6)*1/6)),(ARefTriinv(1,1)*(Coefficients(i,2)+2*Coefficients(i,4)*2/3+Coefficients(i,6)*1/6)+ARefTriinv(2,1)*(Coefficients(i,3)+2*Coefficients(i,5)*1/6+Coefficients(i,6)*2/3))*(ARefTriinv(1,1)*(Coefficients(j,2)+2*Coefficients(j,4)*2/3+Coefficients(j,6)*1/6)+ARefTriinv(2,1)*(Coefficients(j,3)+2*Coefficients(j,5)*1/6+Coefficients(j,6)*2/3))]);
        end
    end
    for i =1:2:5
        for j = 2:2:6
             ALdiv(i,j)=1/6*B*det(ARefTri)*sum([(ARefTriinv(1,1)*(Coefficients(i,2)+2*Coefficients(i,4)*1/6+Coefficients(i,6)*1/6)+ARefTriinv(2,1)*(Coefficients(i,3)+2*Coefficients(i,5)*1/6+Coefficients(i,6)*1/6))*(ARefTriinv(1,2)*(Coefficients(j,2)+2*Coefficients(j,4)*1/6+Coefficients(j,6)*1/6)+ARefTriinv(2,2)*(Coefficients(j,3)+2*Coefficients(j,5)*1/6+Coefficients(j,6)*1/6)),(ARefTriinv(1,1)*(Coefficients(i,2)+2*Coefficients(i,4)*1/6+Coefficients(i,6)*2/3)+ARefTriinv(2,1)*(Coefficients(i,3)+2*Coefficients(i,5)*2/3+Coefficients(i,6)*1/6))*(ARefTriinv(1,2)*(Coefficients(j,2)+2*Coefficients(j,4)*1/6+Coefficients(j,6)*2/3)+ARefTriinv(2,2)*(Coefficients(j,3)+2*Coefficients(j,5)*2/3+Coefficients(j,6)*1/6)),(ARefTriinv(1,1)*(Coefficients(i,2)+2*Coefficients(i,4)*2/3+Coefficients(i,6)*1/6)+ARefTriinv(2,1)*(Coefficients(i,3)+2*Coefficients(i,5)*1/6+Coefficients(i,6)*2/3))*(ARefTriinv(1,2)*(Coefficients(j,2)+2*Coefficients(j,4)*2/3+Coefficients(j,6)*1/6)+ARefTriinv(2,2)*(Coefficients(j,3)+2*Coefficients(j,5)*1/6+Coefficients(j,6)*2/3))]);    
        end
    end
    for i =2:2:6
        for j = 1:2:5
             ALdiv(i,j)=1/6*B*det(ARefTri)*sum([(ARefTriinv(1,2)*(Coefficients(i,2)+2*Coefficients(i,4)*1/6+Coefficients(i,6)*1/6)+ARefTriinv(2,2)*(Coefficients(i,3)+2*Coefficients(i,5)*1/6+Coefficients(i,6)*1/6))*(ARefTriinv(1,1)*(Coefficients(j,2)+2*Coefficients(j,4)*1/6+Coefficients(j,6)*1/6)+ARefTriinv(2,1)*(Coefficients(j,3)+2*Coefficients(j,5)*1/6+Coefficients(j,6)*1/6)),(ARefTriinv(1,2)*(Coefficients(i,2)+2*Coefficients(i,4)*1/6+Coefficients(i,6)*2/3)+ARefTriinv(2,2)*(Coefficients(i,3)+2*Coefficients(i,5)*2/3+Coefficients(i,6)*1/6))*(ARefTriinv(1,1)*(Coefficients(j,2)+2*Coefficients(j,4)*1/6+Coefficients(j,6)*2/3)+ARefTriinv(2,1)*(Coefficients(j,3)+2*Coefficients(j,5)*2/3+Coefficients(j,6)*1/6)),(ARefTriinv(1,2)*(Coefficients(i,2)+2*Coefficients(i,4)*2/3+Coefficients(i,6)*1/6)+ARefTriinv(2,2)*(Coefficients(i,3)+2*Coefficients(i,5)*1/6+Coefficients(i,6)*2/3))*(ARefTriinv(1,1)*(Coefficients(j,2)+2*Coefficients(j,4)*2/3+Coefficients(j,6)*1/6)+ARefTriinv(2,1)*(Coefficients(j,3)+2*Coefficients(j,5)*1/6+Coefficients(j,6)*2/3))]);
        end
    end
    for i =2:2:6
        for j = 2:2:6
             ALdiv(i,j)=1/6*B*det(ARefTri)*sum([(ARefTriinv(1,2)*(Coefficients(i,2)+2*Coefficients(i,4)*1/6+Coefficients(i,6)*1/6)+ARefTriinv(2,2)*(Coefficients(i,3)+2*Coefficients(i,5)*1/6+Coefficients(i,6)*1/6))*(ARefTriinv(1,2)*(Coefficients(j,2)+2*Coefficients(j,4)*1/6+Coefficients(j,6)*1/6)+ARefTriinv(2,2)*(Coefficients(j,3)+2*Coefficients(j,5)*1/6+Coefficients(j,6)*1/6)),(ARefTriinv(1,2)*(Coefficients(i,2)+2*Coefficients(i,4)*1/6+Coefficients(i,6)*2/3)+ARefTriinv(2,2)*(Coefficients(i,3)+2*Coefficients(i,5)*2/3+Coefficients(i,6)*1/6))*(ARefTriinv(1,2)*(Coefficients(j,2)+2*Coefficients(j,4)*1/6+Coefficients(j,6)*2/3)+ARefTriinv(2,2)*(Coefficients(j,3)+2*Coefficients(j,5)*2/3+Coefficients(j,6)*1/6)),(ARefTriinv(1,2)*(Coefficients(i,2)+2*Coefficients(i,4)*2/3+Coefficients(i,6)*1/6)+ARefTriinv(2,2)*(Coefficients(i,3)+2*Coefficients(i,5)*1/6+Coefficients(i,6)*2/3))*(ARefTriinv(1,2)*(Coefficients(j,2)+2*Coefficients(j,4)*2/3+Coefficients(j,6)*1/6)+ARefTriinv(2,2)*(Coefficients(j,3)+2*Coefficients(j,5)*1/6+Coefficients(j,6)*2/3))]);
        end
    end
    
    

    
    ALreg=zeros(6,6);
    dphix = @(x,y,i) Coefficients(i,2)+2*Coefficients(i,4)*x+Coefficients(i,6)*y;
    dphiy = @(x,y,i) Coefficients(i,3)+2*Coefficients(i,5)*y+Coefficients(i,6)*x;
    for i =1:2:5
        for j=1:2:5
            ALreg(i,j)=1/6*det(ARefTri)*sum([(ARefTriinv(1,1)*dphix(1/6,1/6,i)+ARefTriinv(2,1)*dphiy(1/6,1/6,i))*(ARefTriinv(1,1)*dphix(1/6,1/6,j)+ARefTriinv(2,1)*dphiy(1/6,1/6,j))+1/2*(dphiy(1/6,1/6,i)*ARefTriinv(2,2)+dphix(1/6,1/6,i)*ARefTriinv(1,2))*(dphiy(1/6,1/6,j)*ARefTriinv(2,2)+dphix(1/6,1/6,j)*ARefTriinv(1,2)),(ARefTriinv(1,1)*dphix(2/3,1/6,i)+ARefTriinv(2,1)*dphiy(2/3,1/6,i))*(ARefTriinv(1,1)*dphix(2/3,1/6,j)+ARefTriinv(2,1)*dphiy(2/3,1/6,j))+1/2*(dphiy(2/3,1/6,i)*ARefTriinv(2,2)+dphix(2/3,1/6,i)*ARefTriinv(1,2))*(dphiy(2/3,1/6,j)*ARefTriinv(2,2)+dphix(2/3,1/6,j)*ARefTriinv(1,2)),(ARefTriinv(1,1)*dphix(1/6,2/3,i)+ARefTriinv(2,1)*dphiy(1/6,2/3,i))*(ARefTriinv(1,1)*dphix(1/6,2/3,j)+ARefTriinv(2,1)*dphiy(1/6,2/3,j))+1/2*(dphiy(1/6,2/3,i)*ARefTriinv(2,2)+dphix(1/6,2/3,i)*ARefTriinv(1,2))*(dphiy(1/6,2/3,j)*ARefTriinv(2,2)+dphix(1/6,2/3,j)*ARefTriinv(1,2))]);
        end
    end
    for i =2:2:6
        for j=2:2:6
            ALreg(i,j)=1/6*det(ARefTri)*sum([(dphix(1/6,1/6,i)*ARefTriinv(1,2)+dphiy(1/6,1/6,i)*ARefTriinv(2,2))*(dphix(1/6,1/6,j)*ARefTriinv(1,2)+dphiy(1/6,1/6,j)*ARefTriinv(2,2))+1/2*(dphiy(1/6,1/6,i)*ARefTriinv(2,1)+dphix(1/6,1/6,i)*ARefTriinv(1,1))*(dphiy(1/6,1/6,j)*ARefTriinv(2,1)+dphix(1/6,1/6,j)*ARefTriinv(1,1)),(dphix(2/3,1/6,i)*ARefTriinv(1,2)+dphiy(2/3,1/6,i)*ARefTriinv(2,2))*(dphix(2/3,1/6,j)*ARefTriinv(1,2)+dphiy(2/3,1/6,j)*ARefTriinv(2,2))+1/2*(dphiy(2/3,1/6,i)*ARefTriinv(2,1)+dphix(2/3,1/6,i)*ARefTriinv(1,1))*(dphiy(2/3,1/6,j)*ARefTriinv(2,1)+dphix(2/3,1/6,j)*ARefTriinv(1,1)),(dphix(1/6,2/3,i)*ARefTriinv(1,2)+dphiy(1/6,2/3,i)*ARefTriinv(2,2))*(dphix(1/6,2/3,j)*ARefTriinv(1,2)+dphiy(1/6,2/3,j)*ARefTriinv(2,2))+1/2*(dphiy(1/6,2/3,i)*ARefTriinv(2,1)+dphix(1/6,2/3,i)*ARefTriinv(1,1))*(dphiy(1/6,2/3,j)*ARefTriinv(2,1)+dphix(1/6,2/3,j)*ARefTriinv(1,1))]);
        end
    end
    for i =1:2:5
        for j=2:2:6
            ALreg(i,j)=1/12*det(ARefTri)*sum([(dphiy(1/6,1/6,i)*ARefTriinv(2,2)+dphix(1/6,1/6,i)*ARefTriinv(1,2))*(dphiy(1/6,1/6,j)*ARefTriinv(2,1)+dphix(1/6,1/6,j)*ARefTriinv(1,1)),(dphiy(2/3,1/6,i)*ARefTriinv(2,2)+dphix(2/3,1/6,i)*ARefTriinv(1,2))*(dphiy(2/3,1/6,j)*ARefTriinv(2,1)+dphix(2/3,1/6,j)*ARefTriinv(1,1)),(dphiy(1/6,2/3,i)*ARefTriinv(2,2)+dphix(1/6,2/3,i)*ARefTriinv(1,2))*(dphiy(1/6,2/3,j)*ARefTriinv(2,1)+dphix(1/6,2/3,j)*ARefTriinv(1,1))]);
        end
    end
    for i =2:2:6
        for j=1:2:5
            ALreg(i,j)=1/12*det(ARefTri)*sum([(dphiy(1/6,1/6,j)*ARefTriinv(2,2)+dphix(1/6,1/6,j)*ARefTriinv(1,2))*(dphiy(1/6,1/6,i)*ARefTriinv(2,1)+dphix(1/6,1/6,i)*ARefTriinv(1,1)),(dphiy(2/3,1/6,j)*ARefTriinv(2,2)+dphix(2/3,1/6,j)*ARefTriinv(1,2))*(dphiy(2/3,1/6,i)*ARefTriinv(2,1)+dphix(2/3,1/6,i)*ARefTriinv(1,1)),(dphiy(1/6,2/3,j)*ARefTriinv(2,2)+dphix(1/6,2/3,j)*ARefTriinv(1,2))*(dphiy(1/6,2/3,i)*ARefTriinv(2,1)+dphix(1/6,2/3,i)*ARefTriinv(1,1))]);
        end
    end
    
    
    
    A(nodes,nodes) = A(nodes,nodes) + ALdiv+C*ALreg;
    % b
    GaussNodes= [ARefTri*[1/6;1/6]+bRefTri,ARefTri*[2/3;1/6]+bRefTri,ARefTri*[1/6;2/3]+bRefTri];
    fGauss = [f(GaussNodes(1,1),GaussNodes(2,1)),f(GaussNodes(1,2),GaussNodes(2,2)),f(GaussNodes(1,3),GaussNodes(2,3))];
    VectCoeff1=[1,vertices(1,1),vertices(1,2);1,vertices(2,1),vertices(2,2);1,vertices(3,1),vertices(3,2)]\(fVector([nodes(1),nodes(3),nodes(5)]));
    VectCoeff2=[1,vertices(1,1),vertices(1,2);1,vertices(2,1),vertices(2,2);1,vertices(3,1),vertices(3,2)]\(fVector([nodes(2),nodes(4),nodes(6)]));
    fVectGauss1 = [VectCoeff1'*[1;GaussNodes(:,1)],VectCoeff1'*[1;GaussNodes(:,2)],VectCoeff1'*[1;GaussNodes(:,3)]];
    fVectGauss2 = [VectCoeff2'*[1;GaussNodes(:,1)],VectCoeff2'*[1;GaussNodes(:,2)],VectCoeff2'*[1;GaussNodes(:,3)]];
    hVectGauss = [hVect'*[1;GaussNodes(:,1)],hVect'*[1;GaussNodes(:,2)],hVect'*[1;GaussNodes(:,3)]];
    BL=zeros(6,1);
    for j= 1:2:5
    BL(j)=sum((fGauss(1,:)+fVectGauss1).*GaussValues(j,:))+(ARefTriinv(1,1)*Coefficients(j,2)+ARefTriinv(2,1)*Coefficients(j,3))*sum(hVectGauss);
    end
    for j=2:2:6
    BL(j)=sum((fGauss(2,:)+fVectGauss2).*GaussValues(j,:))+(ARefTriinv(1,2)*Coefficients(j,2)+ARefTriinv(2,2)*Coefficients(j,3))*sum(hVectGauss);
    end
    b(nodes) = b(nodes) + 1/6*det(ARefTri)*BL;
    end
    
    
    %Boundary
    for i=1:length(Dirichlet)
        A(Dirichlet(i),:)=0; A(Dirichlet(i),Dirichlet(i))=1;
        b(Dirichlet(i))=DirichletValue;
    end
    
%    for j = 1:size(Neumann,1)
%         nodes = Neumann(j,:);
%         vertices = Coordinates(nodes,:);
%         E = vertices(2,:) - vertices(1,:);
%         b(nodes) = b(nodes) + norm(E)/2 * g(sum(vertices)/2);
%   end
    u=A\b;
    end