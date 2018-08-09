%% Linear Strain equation solver  C(e(u):e(v))+B(nabla*u,nabla*v)=(f,v) +(h,nabla * v)   in 2D domain with 2D Range
function [u,A] = FEMLinStrain2Dvectorial(Coordinates,Elements,C,B,Dirichlet,Neumann,g,f,h,fVector,firstrun,Apre,Coefficients,GaussValues)
    NN=length(Coordinates(:,1));
    b=zeros(2*NN,1);
    number_of_elements=length(Elements(:,1));
    Coefficients = [Coefficients(1,:);Coefficients(1,:);Coefficients(2,:);Coefficients(2,:);Coefficients(3,:);Coefficients(3,:)];
    GaussValues = [GaussValues(1,:);GaussValues(1,:);GaussValues(2,:);GaussValues(2,:);GaussValues(3,:);GaussValues(3,:)];
    if firstrun == 1
    A=sparse(2*NN,2*NN);
   
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
             %ALdiv(i,j) = 1/2*B*[Coefficients(i,2),Coefficients(i,3)]*(ARefTriInv*ARefTriInv')*[Coefficients(j,2);Coefficients(j,3)]*det(ARefTri);
             ALdiv(i,j)=1/2*B*(ARefTriinv(1,1)*Coefficients(i,2)+ARefTriinv(2,1)*Coefficients(i,3))*(ARefTriinv(1,1)*Coefficients(j,2)+ARefTriinv(2,1)*Coefficients(j,3))*det(ARefTri);
        end
    end
    for i =1:2:5
        for j = 2:2:6
             %ALdiv(i,j) = 1/2*B*[Coefficients(i,2),Coefficients(i,3)]*(ARefTriInv*ARefTriInv')*[Coefficients(j,2);Coefficients(j,3)]*det(ARefTri);
             ALdiv(i,j)=1/2*B*(ARefTriinv(1,1)*Coefficients(i,2)+ARefTriinv(2,1)*Coefficients(i,3))*(ARefTriinv(1,2)*Coefficients(j,2)+ARefTriinv(2,2)*Coefficients(j,3))*det(ARefTri);      
        end
    end
    for i =2:2:6
        for j = 1:2:5
             %ALdiv(i,j) = 1/2*B*[Coefficients(i,2),Coefficients(i,3)]*(ARefTriInv*ARefTriInv')*[Coefficients(j,2);Coefficients(j,3)]*det(ARefTri);
             ALdiv(i,j)=1/2*B*(ARefTriinv(1,2)*Coefficients(i,2)+ARefTriinv(2,2)*Coefficients(i,3))*(ARefTriinv(1,1)*Coefficients(j,2)+ARefTriinv(2,1)*Coefficients(j,3))*det(ARefTri);
        end
    end
    for i =2:2:6
        for j = 2:2:6
             %ALdiv(i,j) = 1/2*B*[Coefficients(i,2),Coefficients(i,3)]*(ARefTriInv*ARefTriInv')*[Coefficients(j,2);Coefficients(j,3)]*det(ARefTri);
             ALdiv(i,j)=1/2*B*(ARefTriinv(1,2)*Coefficients(i,2)+ARefTriinv(2,2)*Coefficients(i,3))*(ARefTriinv(1,2)*Coefficients(j,2)+ARefTriinv(2,2)*Coefficients(j,3))*det(ARefTri);
        end
    end
    
    

    
    ALreg=zeros(6,6);
    for i =1:2:5
        for j=1:2:5
            ALreg(i,j)=1/2*det(ARefTri)*((Coefficients(i,2)*ARefTriinv(1,1)+Coefficients(i,3)*ARefTriinv(2,1))*(Coefficients(j,2)*ARefTriinv(1,1)+Coefficients(j,3)*ARefTriinv(2,1))+1/2*(Coefficients(i,3)*ARefTriinv(2,2)+Coefficients(i,2)*ARefTriinv(1,2))*(Coefficients(j,3)*ARefTriinv(2,2)+Coefficients(j,2)*ARefTriinv(1,2)));
        end
    end
    for i =2:2:6
        for j=2:2:6
            ALreg(i,j)=1/2*det(ARefTri)*((Coefficients(i,2)*ARefTriinv(1,2)+Coefficients(i,3)*ARefTriinv(2,2))*(Coefficients(j,2)*ARefTriinv(1,2)+Coefficients(j,3)*ARefTriinv(2,2))+1/2*(Coefficients(i,3)*ARefTriinv(2,1)+Coefficients(i,2)*ARefTriinv(1,1))*(Coefficients(j,3)*ARefTriinv(2,1)+Coefficients(j,2)*ARefTriinv(1,1)));
        end
    end
    for i =1:2:5
        for j=2:2:6
            ALreg(i,j)=1/4*det(ARefTri)*((Coefficients(i,3)*ARefTriinv(2,2)+Coefficients(i,2)*ARefTriinv(1,2))*(Coefficients(j,3)*ARefTriinv(2,1)+Coefficients(j,2)*ARefTriinv(1,1)));
        end
    end
    for i =2:2:6
        for j=1:2:5
            ALreg(i,j)=1/4*det(ARefTri)*((Coefficients(j,3)*ARefTriinv(2,2)+Coefficients(j,2)*ARefTriinv(1,2))*(Coefficients(i,3)*ARefTriinv(2,1)+Coefficients(i,2)*ARefTriinv(1,1)));
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
                A(Dirichlet(i,1),:)=0; A(Dirichlet(i,1),Dirichlet(i,1))=1;
                b(Dirichlet(i,1))=Dirichlet(i,2);
    end
    end
    
    if firstrun ~= 1
    for k = 1:number_of_elements
    nodes = Elements(k,:);
    vertices = Coordinates(nodes,:);
    hVect=[1,vertices(1,1),vertices(1,2);1,vertices(2,1),vertices(2,2);1,vertices(3,1),vertices(3,2)]\h(nodes);
    nodes=[2*nodes(1)-1,2*nodes(1),2*nodes(2)-1,2*nodes(2),2*nodes(3)-1,2*nodes(3)];
    
    
    %Transformation to refTriangle
    [ARefTri,bRefTri]=RefTriangleMapInv(vertices(:,1),vertices(:,2));
    [ARefTriinv,~]=RefTriangleMap(vertices(:,1),vertices(:,2));
    %ARefTriinv=ARefTri^(-1);
        
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
    for i=1:length(Dirichlet)       
           b(Dirichlet(i,1))=Dirichlet(i,2);
    end  
    A=Apre;
    end
    u=A\b;
    end