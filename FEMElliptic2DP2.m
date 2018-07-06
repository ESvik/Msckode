%% Elliptic equation solver  B(nabla*u,nabla*v)=(f,v)   in 2D domain using P2 elements.
function [u,A] = FEMElliptic2DP2(Coordinates,Elements,B,Dirichlet,f,Coefficients,GaussValues,NN)
    b=zeros((2*sqrt(NN)-1)^2,1);
    number_of_elements=length(Elements(:,1));
    
    %Assembly
        A=sparse((2*sqrt(NN)-1)^2,(2*sqrt(NN)-1)^2);
        for k = 1:number_of_elements
            nodes = Elements(k,:);
            originalnodes=[nodes(1),nodes(2),nodes(3)];
            vertices = Coordinates(originalnodes,:);    
    
            %Transformation to refTriangle
            [ARefTri,bRefTri]=RefTriangleMapInv(vertices(:,1),vertices(:,2));
            [ARefTriinv,~]=RefTriangleMap(vertices(:,1),vertices(:,2));
            
            dphix = @(x,y,i) Coefficients(i,2)+2*Coefficients(i,4)*x+Coefficients(i,6)*y;
            dphiy = @(x,y,i) Coefficients(i,3)+2*Coefficients(i,5)*y+Coefficients(i,6)*x;
            ALdivpre = @(x,y,i,j) (ARefTriinv(1,1)*dphix(x,y,i) + ARefTriinv(2,1)*dphiy(x,y,i))*(ARefTriinv(1,1)*dphix(x,y,j) + ARefTriinv(2,1)*dphiy(x,y,j))+(ARefTriinv(1,2)*dphix(x,y,i) + ARefTriinv(2,2)*dphiy(x,y,i))*(ARefTriinv(1,2)*dphix(x,y,j) + ARefTriinv(2,2)*dphiy(x,y,j));
            %A
            ALdiv=zeros(6,6);
            for i =1:6
                for j = 1:6
                    ALdiv(i,j)=B*1/6*(ALdivpre(1/6,1/6,i,j)+ALdivpre(1/6,2/3,i,j)+ALdivpre(2/3,1/6,i,j))*det(ARefTri);
                end
            end
            
          
            A(nodes,nodes) = A(nodes,nodes) + ALdiv;
            % b
            GaussNodes= [ARefTri*[1/6;1/6]+bRefTri,ARefTri*[2/3;1/6]+bRefTri,ARefTri*[1/6;2/3]+bRefTri];
            fGauss = [f(GaussNodes(1,1),GaussNodes(2,1)),f(GaussNodes(1,2),GaussNodes(2,2)),f(GaussNodes(1,3),GaussNodes(2,3))];
            BL=zeros(6,1);
            for j= 1:6
                BL(j)=sum((fGauss(1,:)).*GaussValues(j,:));
            end
            b(nodes) = b(nodes) + 1/6*det(ARefTri)*BL;
    
            %Boundary
            for i=1:length(Dirichlet(:,1))
                A(Dirichlet(i,1),:)=0; A(Dirichlet(i,1),Dirichlet(i,1))=1;
                b(Dirichlet(i,1))=Dirichlet(i,2);
            end
    
        end
    
   
    u=A\b;
    
    end