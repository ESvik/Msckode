%FEMElliptic2DP1

% Solver for B(u',v')=(f,v) in 2D domain

function u = FEMElliptic2DP1(Coordinates,Elements,B,Dirichlet,DirichletValue,f,Coefficients,GaussValues)
    NN=length(Coordinates(:,1));
    A=sparse(NN,NN);
    b=zeros(NN,1);
    number_of_elements=length(Elements(:,1));
    
    %Assembly
   
    
        for k = 1:number_of_elements
            nodes = Elements(k,:);
            vertices = Coordinates(nodes,:);
     
            %Transformation to refTriangle
            [ARefTri,bRefTri]=RefTriangleMapInv(vertices(:,1),vertices(:,2));
            ARefTriInv = ARefTri^(-1);
            ALdiv=zeros(3,3);
            for i =1:3
                for j = 1:3
                    ALdiv(i,j) = 1/2*B*[Coefficients(i,2),Coefficients(i,3)]*(ARefTriInv*ARefTriInv')*[Coefficients(j,2);Coefficients(j,3)]*det(ARefTri);
                end
            end
               A(nodes,nodes)=A(nodes,nodes)+ALdiv;
          % b
              GaussNodes= [ARefTri*[1/6;1/6]+bRefTri,ARefTri*[2/3;1/6]+bRefTri,ARefTri*[1/6;2/3]+bRefTri];
             
            fGauss = [f(GaussNodes(1,1),GaussNodes(2,1)),f(GaussNodes(1,2),GaussNodes(2,2)),f(GaussNodes(1,3),GaussNodes(2,3))];
            BL=zeros(3,1);
            for j= 1:3
                BL(j) = det(ARefTri)*sum((fGauss).*GaussValues(j,:));
            end
            b(nodes) = b(nodes) + BL;
    
    
    
            %Boundary
            for i=1:length(Dirichlet)
                A(Dirichlet(i),:)=0; A(Dirichlet(i),Dirichlet(i))=1;
                b(Dirichlet(i))=DirichletValue;
            end
    
          
        end
        %[L,U]=lu(A);
   
    
   
    u=A\b;
        
    end