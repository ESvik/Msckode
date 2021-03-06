% Solver for C(u,v)+B(u',v')=(f,v) in 2D domain

function [u,A] = FEMParabolic2DP1Mscheme(Coordinates,Elements,kappa,tau,Dirichlet,DirichletValue,Neumann,g,f,fVector,firstrun,Apre,Coefficients,GaussValues,bprime,u_prev)
    NN=length(Coordinates(:,1));
    A=sparse(NN,NN);
    b=zeros(NN,1);
    number_of_elements=length(Elements(:,1));
%     %Basis and its evaluation points
%     Coefficients= zeros(3,3);
%     MatOfCoord=[1,0,0;1,1,0;1,0,1];
%     for i =1:3
%         bb=zeros(3,1);
%         bb(i)=1;
%         Coefficients(i,:)=MatOfCoord\bb;
%     end
%     GaussValues = zeros(3,3);
%     for i = 1:3
%         GaussValues(i,1) = 1/6*(Coefficients(i,1) + Coefficients(i,2)*1/6 + Coefficients(i,3)*1/6);
%         GaussValues(i,2) = 1/6*(Coefficients(i,1) + Coefficients(i,2)*2/3 + Coefficients(i,3)*1/6);
%         GaussValues(i,3) = 1/6*(Coefficients(i,1) + Coefficients(i,2)*1/6 + Coefficients(i,3)*2/3);
%     end
    
    %Assembly
   
    
    %A
    if firstrun == 1
        for k = 1:number_of_elements
            nodes = Elements(k,:);
            vertices = Coordinates(nodes,:);
            Lb=max(bprime(u_prev(nodes)));
            bm=min(bprime(u_prev(nodes)));
            gamma = max(0,1+(Lb*bm-sqrt(Lb^2*bm^2+4*Lb*bm^2*kappa*tau+4*Lb*bm*kappa^2*tau^2))/(4*bm*kappa*tau));
            L=Lb/(2*(1-gamma));
     
            %Transformation to refTriangle
            [ARefTri,bRefTri]=RefTriangleMapInv(vertices(:,1),vertices(:,2));
            ARefTriInv = ARefTri^(-1);
            ALdiv=zeros(3,3);
            for i =1:3
                for j = 1:3
                    ALdiv(i,j) = 1/2*kappa*tau*[Coefficients(i,2),Coefficients(i,3)]*(ARefTriInv*ARefTriInv')*[Coefficients(j,2);Coefficients(j,3)]*det(ARefTri);
                end
            end
    
            ALreg=zeros(3,3);
            for i =1:3
                for j=1:3
                    ALreg(i,j)=1/6*L*det(ARefTri)*sum(GaussValues(i,:).*GaussValues(j,:));
                end
            end
            A(nodes,nodes) = A(nodes,nodes) + ALdiv+ALreg;
    
            % b
            GaussNodes= [ARefTri*[1/6;1/6]+bRefTri,ARefTri*[2/3;1/6]+bRefTri,ARefTri*[1/6;2/3]+bRefTri];
  
            fGauss = [f(GaussNodes(1,1),GaussNodes(2,1)),f(GaussNodes(1,2),GaussNodes(2,2)),f(GaussNodes(1,3),GaussNodes(2,3))];
            VectCoeff=[1,vertices(1,1),vertices(1,2);1,vertices(2,1),vertices(2,2 );1,vertices(3,1),vertices(3,2)]\(fVector(nodes));
            uprevCoeff=[1,vertices(1,1),vertices(1,2);1,vertices(2,1),vertices(2,2 );1,vertices(3,1),vertices(3,2)]\(u_prev(nodes));
            fVectGauss = [VectCoeff'*[1;GaussNodes(:,1)],VectCoeff'*[1;GaussNodes(:,2)],VectCoeff'*[1;GaussNodes(:,3)]];
            uVectGauss = [uprevCoeff'*[1;GaussNodes(:,1)],uprevCoeff'*[1;GaussNodes(:,2)],uprevCoeff'*[1;GaussNodes(:,3)]];
            BL=zeros(3,1);
            for j= 1:3
                BL(j) = 1/6*det(ARefTri)*sum((fGauss+fVectGauss+L*uVectGauss).*GaussValues(j,:));
            end
            b(nodes) = b(nodes) + BL;
    
    
    
            %Boundary
            for i=1:length(Dirichlet)
                A(Dirichlet(i),:)=0; A(Dirichlet(i),Dirichlet(i))=1;
                b(Dirichlet(i))=DirichletValue;
            end
    
            for j = 1:size(Neumann,1)
                nodes = Neumann(j,:);
                vertices = Coordinates(nodes,:);
                E = vertices(2,:) - vertices(1,:);
                b(nodes) = b(nodes) + norm(E)/2 * g(sum(vertices)/2);
            end
        end
        %[L,U]=lu(A);
    end
    
    
    if firstrun ~= 1
        for k = 1:number_of_elements
            nodes = Elements(k,:);
            vertices = Coordinates(nodes,:);
            Lb=max(bprime(u_prev(nodes)));
            bm=min(bprime(u_prev(nodes)));
            gamma = max(0,1+(Lb*bm-sqrt(Lb^2*bm^2+4*Lb*bm^2*kappa*tau+4*Lb*bm*kappa^2*tau^2))/(4*bm*kappa*tau));
            L=Lb/(2*(1-gamma));
            L/Lb;
            
            [ARefTri,bRefTri]=RefTriangleMapInv(vertices(:,1),vertices(:,2));
            ARefTriInv = ARefTri^(-1);
            
     
            %Transformation to refTriangle
            GaussNodes= [ARefTri*[1/6;1/6]+bRefTri,ARefTri*[2/3;1/6]+bRefTri,ARefTri*[1/6;2/3]+bRefTri];
  
            fGauss = [f(GaussNodes(1,1),GaussNodes(2,1)),f(GaussNodes(1,2),GaussNodes(2,2)),f(GaussNodes(1,3),GaussNodes(2,3))];
            VectCoeff=[1,vertices(1,1),vertices(1,2);1,vertices(2,1),vertices(2,2 );1,vertices(3,1),vertices(3,2)]\(fVector(nodes));
            uprevCoeff=[1,vertices(1,1),vertices(1,2);1,vertices(2,1),vertices(2,2 );1,vertices(3,1),vertices(3,2)]\(u_prev(nodes));
            fVectGauss = [VectCoeff'*[1;GaussNodes(:,1)],VectCoeff'*[1;GaussNodes(:,2)],VectCoeff'*[1;GaussNodes(:,3)]];
            uVectGauss = [uprevCoeff'*[1;GaussNodes(:,1)],uprevCoeff'*[1;GaussNodes(:,2)],uprevCoeff'*[1;GaussNodes(:,3)]];
            BL=zeros(3,1);
            for j= 1:3
                BL(j) = 1/6*det(ARefTri)*sum((fGauss+fVectGauss+L*uVectGauss).*GaussValues(j,:));
            end
            b(nodes) = b(nodes) + BL;
    
        end
        A=Apre;
        for i=1:length(Dirichlet)
            b(Dirichlet(i))=DirichletValue;
        end
    end
    
    u=A\b;
        
    end