% Solver for C(u,v)+B(nabla*u,nabla*v)=(f,v) +(h,nabla * v)   in 2D domain with 2D Range

function u = FEMParabolic2Dvectorial(Coordinates,Elements,C,B,Dirichlet,DirichletValue,Neumann,g,f,h,fVector)
    NN=length(Coordinates(:,1));
    A=sparse(2*NN,2*NN);
    b=zeros(2*NN,1);
    number_of_elements=length(Elements(:,1));
    
    %Basis and its evaluation points
    Coefficients= zeros(3,3);
    MatOfCoord=[1,0,0;1,1,0;1,0,1];
    for i =1:3
        bb=zeros(3,1);
        bb(i)=1;
        Coefficients(i,:)=MatOfCoord\bb;
    end
    Coefficients=[Coefficients(1,:);Coefficients(1,:);Coefficients(2,:);Coefficients(2,:);Coefficients(3,:);Coefficients(3,:)];
    GaussValues = zeros(3,3);
    for i = 1:3
        GaussValues(i,1) = (Coefficients(i,1) + Coefficients(i,2)*1/6 + Coefficients(i,3)*1/6);
        GaussValues(i,2) = (Coefficients(i,1) + Coefficients(i,2)*2/3 + Coefficients(i,3)*1/6);
        GaussValues(i,3) = (Coefficients(i,1) + Coefficients(i,2)*1/6 + Coefficients(i,3)*2/3);
    end
    GaussValues=[GaussValues(1,:);GaussValues(1,:);GaussValues(2,:);GaussValues(2,:);GaussValues(3,:);GaussValues(3,:)];
    
    %Assembly
    for k = 1:number_of_elements
    nodes = Elements(k,:);
    vertices = Coordinates(nodes,:);
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
            ALreg(i,j)=1/6*C*det(ARefTri)*sum(GaussValues(i,:).*GaussValues(j,:));
        end
    end
     for i =2:2:6
        for j=2:2:6
            ALreg(i,j)=1/6*C*det(ARefTri)*sum(GaussValues(i,:).*GaussValues(j,:));
        end
     end
    
    
    A(nodes,nodes) = A(nodes,nodes) + ALdiv+ALreg;
    % b
    GaussNodes= [ARefTri*[1/6;1/6]+bRefTri,ARefTri*[2/3;1/6]+bRefTri,ARefTri*[1/6;2/3]+bRefTri];
    fGauss = [f(GaussNodes(1,1),GaussNodes(2,1)),f(GaussNodes(1,2),GaussNodes(2,2)),f(GaussNodes(1,3),GaussNodes(2,3))];
    
%     VectCoeff=[1,vertices(1,1),vertices(1,2);1,vertices(2,1),vertices(2,2 );1,vertices(3,1),vertices(3,2)]\(fVector(nodes));
%     fVectGauss = [VectCoeff'*[1;GaussNodes(:,1)],VectCoeff'*[1;GaussNodes(:,2)],VectCoeff'*[1;GaussNodes(:,3)]];
    BL=zeros(6,1);
    for j= 1:2:5
    BL(j)=sum(fGauss(1,:).*GaussValues(j,:));
    end
    for j=2:2:6
    BL(j) =sum(fGauss(2,:).*GaussValues(j,:));
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
    etha=A\b;
    u=zeros(NN,2);
    u(:,1)=etha(1:2:2*NN-1); u(:,2)=etha(2:2:2*NN);
    end