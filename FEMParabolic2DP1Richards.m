% Solver for C(u,v)+B(u',v')=(f,v) in 2D domain
% 1 = L-scheme, 2 = Modified L-scheme, 3 = Newton, , 4 = local L-scheme
function [u,A] = FEMParabolic2DP1Richards(Coordinates,Elements,C,B,Dirichlet,Neumann,g,f,fVector,firstrun,Apre,Coefficients,GaussValues,bprime,u_prev,scheme,modifiedconstant,Kirchhoff,kappa,kappaprime,tau)
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
   
    if Kirchhoff ==1
    %% L-scheme
    if scheme == 1
        fVector = fVector + C*u_prev;
    if firstrun == 1
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
    
            ALreg=zeros(3,3);
            for i =1:3
                for j=1:3
                    ALreg(i,j)=1/6*C*det(ARefTri)*sum(GaussValues(i,:).*GaussValues(j,:));
                end
            end
            A(nodes,nodes) = A(nodes,nodes) + ALdiv+ALreg;
    
            % b
            %biotdisplacement
            %modes=[2*nodes(1)-1,2*nodes(1),2*nodes(2)-1,2*nodes(2),2*nodes(3)-1,2*nodes(3)];
            GaussNodes= [ARefTri*[1/6;1/6]+bRefTri,ARefTri*[2/3;1/6]+bRefTri,ARefTri*[1/6;2/3]+bRefTri];
            
            fGauss = [f(GaussNodes(1,1),GaussNodes(2,1)),f(GaussNodes(1,2),GaussNodes(2,2)),f(GaussNodes(1,3),GaussNodes(2,3))];
            VectCoeff=[1,vertices(1,1),vertices(1,2);1,vertices(2,1),vertices(2,2 );1,vertices(3,1),vertices(3,2)]\(fVector(nodes));
            fVectGauss = [VectCoeff'*[1;GaussNodes(:,1)],VectCoeff'*[1;GaussNodes(:,2)],VectCoeff'*[1;GaussNodes(:,3)]];
            BL=zeros(3,1);
            for j= 1:3
                BL(j) = 1/6*det(ARefTri)*sum((fGauss+fVectGauss).*GaussValues(j,:));
            end
            b(nodes) = b(nodes) + BL;
    
    
    
            %Boundary
            for i=1:length(Dirichlet(:,1))
                A(Dirichlet(i,1),:)=0; A(Dirichlet(i,1),Dirichlet(i,1))=1;
                b(Dirichlet(i,1))=Dirichlet(i,2);
            end
            end
    
            for j = 1:size(Neumann,1)
                nodes = Neumann(j,:);
                vertices = Coordinates(nodes,:);
                E = vertices(2,:) - vertices(1,:);
                b(nodes) = b(nodes) + norm(E)/2 * g(sum(vertices)/2);
            end
        end
        %[L,U]=lu(A);
    
    
    
    if firstrun ~= 1
        for k = 1:number_of_elements
            nodes = Elements(k,:);
            vertices = Coordinates(nodes,:);
            
     
            %Transformation to refTriangle
            [ARefTri,bRefTri]=RefTriangleMapInv(vertices(:,1),vertices(:,2));
            % b
            %biotdisplacement
            GaussNodes= [ARefTri*[1/6;1/6]+bRefTri,ARefTri*[2/3;1/6]+bRefTri,ARefTri*[1/6;2/3]+bRefTri];
            fGauss = [f(GaussNodes(1,1),GaussNodes(2,1)),f(GaussNodes(1,2),GaussNodes(2,2)),f(GaussNodes(1,3),GaussNodes(2,3))];
            VectCoeff=[1,vertices(1,1),vertices(1,2);1,vertices(2,1),vertices(2,2 );1,vertices(3,1),vertices(3,2)]\(fVector(nodes));
            fVectGauss = [VectCoeff'*[1;GaussNodes(:,1)],VectCoeff'*[1;GaussNodes(:,2)],VectCoeff'*[1;GaussNodes(:,3)]];
            BL=zeros(3,1);
            for j= 1:3
                BL(j) = 1/6*det(ARefTri)*sum((fGauss+fVectGauss).*GaussValues(j,:));
            end
            b(nodes) = b(nodes) + BL;
    
       end
        A=Apre;
        for i=1:length(Dirichlet(:,1))
            b(Dirichlet(i,1))=Dirichlet(i,2);
        end
    end
    
    u=A\b;
        
      
    end
    
    %% Modified L-scheme
    if scheme == 2
        for k = 1:number_of_elements
            nodes = Elements(k,:);
            vertices = Coordinates(nodes,:);
    
            m=modifiedconstant;
            %Transformation to refTriangle
            [ARefTri,bRefTri]=RefTriangleMapInv(vertices(:,1),vertices(:,2));
            ARefTriInv = ARefTri^(-1);
            GaussNodes= [ARefTri*[1/6;1/6]+bRefTri,ARefTri*[2/3;1/6]+bRefTri,ARefTri*[1/6;2/3]+bRefTri];
            MM=[max(bprime(u_prev(nodes(1)))+m,2*m),max(bprime(u_prev(nodes(2)))+m,2*m),max(bprime(u_prev(nodes(3)))+m,2*m)];
            MVectCoeff=[1,vertices(1,1),vertices(1,2);1,vertices(2,1),vertices(2,2 );1,vertices(3,1),vertices(3,2)]\(MM');
            M = [MVectCoeff'*[1;GaussNodes(:,1)],MVectCoeff'*[1;GaussNodes(:,2)],MVectCoeff'*[1;GaussNodes(:,3)]];
            
            ALdiv=zeros(3,3);
            for i =1:3
                for j = 1:3
                    ALdiv(i,j) = 1/2*B*[Coefficients(i,2),Coefficients(i,3)]*(ARefTriInv*ARefTriInv')*[Coefficients(j,2);Coefficients(j,3)]*det(ARefTri);
                end
            end
    
            ALreg=zeros(3,3);
            for i =1:3
                for j=1:3
                    ALreg(i,j)=1/6*det(ARefTri)*sum(M.*GaussValues(i,:).*GaussValues(j,:));
                end
            end
            A(nodes,nodes) = A(nodes,nodes) + ALdiv+ALreg;
    
            % b
            %biotdisplacement
            %modes=[2*nodes(1)-1,2*nodes(1),2*nodes(2)-1,2*nodes(2),2*nodes(3)-1,2*nodes(3)];
            
            fGauss = [f(GaussNodes(1,1),GaussNodes(2,1)),f(GaussNodes(1,2),GaussNodes(2,2)),f(GaussNodes(1,3),GaussNodes(2,3))];
            VectCoeff=[1,vertices(1,1),vertices(1,2);1,vertices(2,1),vertices(2,2 );1,vertices(3,1),vertices(3,2)]\(fVector(nodes));
            fVectGauss = [VectCoeff'*[1;GaussNodes(:,1)],VectCoeff'*[1;GaussNodes(:,2)],VectCoeff'*[1;GaussNodes(:,3)]];
            uprevCoeff=[1,vertices(1,1),vertices(1,2);1,vertices(2,1),vertices(2,2 );1,vertices(3,1),vertices(3,2)]\(u_prev(nodes));
            uVectGauss = [uprevCoeff'*[1;GaussNodes(:,1)],uprevCoeff'*[1;GaussNodes(:,2)],uprevCoeff'*[1;GaussNodes(:,3)]];
           
            BL=zeros(3,1);
            for j= 1:3
                BL(j) = 1/6*det(ARefTri)*sum((fGauss+fVectGauss+M.*uVectGauss).*GaussValues(j,:));
            end
            b(nodes) = b(nodes) + BL;
    
    
    
            %Boundary
            for i=1:length(Dirichlet(:,1))
                A(Dirichlet(i,1),:)=0; A(Dirichlet(i,1),Dirichlet(i,1))=1;
                b(Dirichlet(i,1))=Dirichlet(i,2);
            end
    
            for j = 1:size(Neumann,1)
                nodes = Neumann(j,:);
                vertices = Coordinates(nodes,:);
                E = vertices(2,:) - vertices(1,:);
                b(nodes) = b(nodes) + norm(E)/2 * g(sum(vertices)/2);
            end
        end
        %[L,U]=lu(A);
    
    
    
    u=A\b;
    end
    

    %% Newton
    if scheme == 3
        for k = 1:number_of_elements
            nodes = Elements(k,:);
            vertices = Coordinates(nodes,:);
            
    
     
            %Transformation to refTriangle
            [ARefTri,bRefTri]=RefTriangleMapInv(vertices(:,1),vertices(:,2));
            ARefTriInv = ARefTri^(-1);
            GaussNodes= [ARefTri*[1/6;1/6]+bRefTri,ARefTri*[2/3;1/6]+bRefTri,ARefTri*[1/6;2/3]+bRefTri];
            MVectCoeff=[1,vertices(1,1),vertices(1,2);1,vertices(2,1),vertices(2,2 );1,vertices(3,1),vertices(3,2)]\(bprime(u_prev(nodes)));
            M = [MVectCoeff'*[1;GaussNodes(:,1)],MVectCoeff'*[1;GaussNodes(:,2)],MVectCoeff'*[1;GaussNodes(:,3)]];
            
            ALdiv=zeros(3,3);
            for i =1:3
                for j = 1:3
                    ALdiv(i,j) = 1/2*B*[Coefficients(i,2),Coefficients(i,3)]*(ARefTriInv*ARefTriInv')*[Coefficients(j,2);Coefficients(j,3)]*det(ARefTri);
                end
            end
    
            ALreg=zeros(3,3);
            for i =1:3
                for j=1:3
                    ALreg(i,j)=1/6*det(ARefTri)*sum(M.*GaussValues(i,:).*GaussValues(j,:));
                end
            end
            A(nodes,nodes) = A(nodes,nodes) + ALdiv+ALreg;
    
            % b
            %biotdisplacement
            %modes=[2*nodes(1)-1,2*nodes(1),2*nodes(2)-1,2*nodes(2),2*nodes(3)-1,2*nodes(3)];
            
            fGauss = [f(GaussNodes(1,1),GaussNodes(2,1)),f(GaussNodes(1,2),GaussNodes(2,2)),f(GaussNodes(1,3),GaussNodes(2,3))];
            VectCoeff=[1,vertices(1,1),vertices(1,2);1,vertices(2,1),vertices(2,2 );1,vertices(3,1),vertices(3,2)]\(fVector(nodes));
            fVectGauss = [VectCoeff'*[1;GaussNodes(:,1)],VectCoeff'*[1;GaussNodes(:,2)],VectCoeff'*[1;GaussNodes(:,3)]];
            uprevCoeff1=[1,vertices(1,1),vertices(1,2);1,vertices(2,1),vertices(2,2 );1,vertices(3,1),vertices(3,2)]\(u_prev(2*nodes-1));
            uVectGauss1 = [uprevCoeff1'*[1;GaussNodes(:,1)],uprevCoeff1'*[1;GaussNodes(:,2)],uprevCoeff1'*[1;GaussNodes(:,3)]];
            uprevCoeff2=[1,vertices(1,1),vertices(1,2);1,vertices(2,1),vertices(2,2 );1,vertices(3,1),vertices(3,2)]\(u_prev(2*nodes));
            uVectGauss2 = [uprevCoeff2'*[1;GaussNodes(:,1)],uprevCoeff2'*[1;GaussNodes(:,2)],uprevCoeff2'*[1;GaussNodes(:,3)]];
            eps = uVectCoeff1(2)^2+0.5*(uVectCoeff1(3)+uVectCoeff2(2))^2+uVectCoeff2(3)^2;
            div = m
            
            
            
            
            BL=zeros(3,1);
            for j= 1:3
                BL(j) = 1/6*det(ARefTri)*sum((fGauss+fVectGauss+M.*uVectGauss).*GaussValues(j,:));
            end
            b(nodes) = b(nodes) + BL;
    
    
    
            %Boundary
            for i=1:length(Dirichlet(:,1))
                A(Dirichlet(i,1),:)=0; A(Dirichlet(i,1),Dirichlet(i,1))=1;
                b(Dirichlet(i,1))=Dirichlet(i,2);
            end
    
            for j = 1:size(Neumann,1)
                nodes = Neumann(j,:);
                vertices = Coordinates(nodes,:);
                E = vertices(2,:) - vertices(1,:);
                b(nodes) = b(nodes) + norm(E)/2 * g(sum(vertices)/2);
            end
        end
        %[L,U]=lu(A);
    
    

    
    u=A\b;
        
    end

    %% Local L-scheme
    if scheme == 4
        
        for k = 1:number_of_elements
            nodes = Elements(k,:);
            vertices = Coordinates(nodes,:);
            Lb=max(bprime(u_prev(nodes)));
            bm=min(bprime(u_prev(nodes)));
            gamma = max(0,1+(Lb*bm-sqrt(Lb^2*bm^2+4*Lb*bm^2*B+4*Lb*bm*B^2))/(4*bm*B));
            L=Lb/(2*(1-gamma));
     
            %Transformation to refTriangle
            [ARefTri,bRefTri]=RefTriangleMapInv(vertices(:,1),vertices(:,2));
            ARefTriInv = ARefTri^(-1);
            ALdiv=zeros(3,3);
            for i =1:3
                for j = 1:3
                    ALdiv(i,j) = 1/2*B*[Coefficients(i,2),Coefficients(i,3)]*(ARefTriInv*ARefTriInv')*[Coefficients(j,2);Coefficients(j,3)]*det(ARefTri);
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
            for i=1:length(Dirichlet(:,1))
                A(Dirichlet(i,1),:)=0; A(Dirichlet(i,1),Dirichlet(i,1))=1;
                b(Dirichlet(i,1))=Dirichlet(i,2);
            end
    
            for j = 1:size(Neumann,1)
                nodes = Neumann(j,:);
                vertices = Coordinates(nodes,:);
                E = vertices(2,:) - vertices(1,:);
                b(nodes) = b(nodes) + norm(E)/2 * g(sum(vertices)/2);
            end
        end
        %[L,U]=lu(A);
    
    
    
   
    u=A\b;
    end
    end

%% -------------------------------------------------------------------------------
    
    if Kirchhoff ==2
    %% L-scheme
    if scheme == 1
        fVector = fVector + C*u_prev;
        for k = 1:number_of_elements
            nodes = Elements(k,:);
            vertices = Coordinates(nodes,:);
    
            [ARefTri,bRefTri]=RefTriangleMapInv(vertices(:,1),vertices(:,2));
            ARefTriInv = ARefTri^(-1);
            GaussNodes= [ARefTri*[1/6;1/6]+bRefTri,ARefTri*[2/3;1/6]+bRefTri,ARefTri*[1/6;2/3]+bRefTri];
            
            KappaVectCoeff=[1,vertices(1,1),vertices(1,2);1,vertices(2,1),vertices(2,2 );1,vertices(3,1),vertices(3,2)]\(kappa(u_prev(nodes)));
            Kappa = [KappaVectCoeff'*[1;GaussNodes(:,1)],KappaVectCoeff'*[1;GaussNodes(:,2)],KappaVectCoeff'*[1;GaussNodes(:,3)]];

            ALdiv=zeros(3,3);
            for i =1:3
                for j = 1:3
                    ALdiv(i,j) = 1/6*tau*[Coefficients(i,2),Coefficients(i,3)]*(ARefTriInv*ARefTriInv')*[Coefficients(j,2);Coefficients(j,3)]*det(ARefTri)*sum(Kappa);
                end
            end
    
            ALreg=zeros(3,3);
            for i =1:3
                for j=1:3
                    ALreg(i,j)=1/6*C*det(ARefTri)*sum(GaussValues(i,:).*GaussValues(j,:));
                end
            end
            A(nodes,nodes) = A(nodes,nodes) + ALdiv+ALreg;
    
            % b
            %biotdisplacement
            %modes=[2*nodes(1)-1,2*nodes(1),2*nodes(2)-1,2*nodes(2),2*nodes(3)-1,2*nodes(3)];
            GaussNodes= [ARefTri*[1/6;1/6]+bRefTri,ARefTri*[2/3;1/6]+bRefTri,ARefTri*[1/6;2/3]+bRefTri];
            
            fGauss = [f(GaussNodes(1,1),GaussNodes(2,1)),f(GaussNodes(1,2),GaussNodes(2,2)),f(GaussNodes(1,3),GaussNodes(2,3))];
            VectCoeff=[1,vertices(1,1),vertices(1,2);1,vertices(2,1),vertices(2,2 );1,vertices(3,1),vertices(3,2)]\(fVector(nodes));
            fVectGauss = [VectCoeff'*[1;GaussNodes(:,1)],VectCoeff'*[1;GaussNodes(:,2)],VectCoeff'*[1;GaussNodes(:,3)]];
            BL=zeros(3,1);
            for j= 1:3
                BL(j) = 1/6*det(ARefTri)*sum((fGauss+fVectGauss).*GaussValues(j,:));
            end
            b(nodes) = b(nodes) + BL;
    
    
    
            %Boundary
            for i=1:length(Dirichlet(:,1))
                A(Dirichlet(i,1),:)=0; A(Dirichlet(i,1),Dirichlet(i,1))=1;
                b(Dirichlet(i,1))=Dirichlet(i,2);
            end
    
            for j = 1:size(Neumann,1)
                nodes = Neumann(j,:);
                vertices = Coordinates(nodes,:);
                E = vertices(2,:) - vertices(1,:);
                b(nodes) = b(nodes) + norm(E)/2 * g(sum(vertices)/2);
            end
        end
        %[L,U]=lu(A);
   
    u=A\b;
        
    end
    
    %% Modified L-scheme
    if scheme == 2
        for k = 1:number_of_elements
            nodes = Elements(k,:);
            vertices = Coordinates(nodes,:);
    
            m=modifiedconstant;
            %Transformation to refTriangle
            [ARefTri,bRefTri]=RefTriangleMapInv(vertices(:,1),vertices(:,2));
            ARefTriInv = ARefTri^(-1);
            GaussNodes= [ARefTri*[1/6;1/6]+bRefTri,ARefTri*[2/3;1/6]+bRefTri,ARefTri*[1/6;2/3]+bRefTri];
            MM=[max(bprime(u_prev(nodes(1)))+m,2*m),max(bprime(u_prev(nodes(2)))+m,2*m),max(bprime(u_prev(nodes(3)))+m,2*m)];
            MVectCoeff=[1,vertices(1,1),vertices(1,2);1,vertices(2,1),vertices(2,2 );1,vertices(3,1),vertices(3,2)]\(MM');
            M = [MVectCoeff'*[1;GaussNodes(:,1)],MVectCoeff'*[1;GaussNodes(:,2)],MVectCoeff'*[1;GaussNodes(:,3)]];
            
            KappaVectCoeff=[1,vertices(1,1),vertices(1,2);1,vertices(2,1),vertices(2,2 );1,vertices(3,1),vertices(3,2)]\(kappa(u_prev(nodes)));
            Kappa = [KappaVectCoeff'*[1;GaussNodes(:,1)],KappaVectCoeff'*[1;GaussNodes(:,2)],KappaVectCoeff'*[1;GaussNodes(:,3)]];

            ALdiv=zeros(3,3);
            for i =1:3
                for j = 1:3
                    ALdiv(i,j) = 1/6*tau*[Coefficients(i,2),Coefficients(i,3)]*(ARefTriInv*ARefTriInv')*[Coefficients(j,2);Coefficients(j,3)]*det(ARefTri)*sum(Kappa);
                end
            end
    
            ALreg=zeros(3,3);
            for i =1:3
                for j=1:3
                    ALreg(i,j)=1/6*det(ARefTri)*sum(M.*GaussValues(i,:).*GaussValues(j,:));
                end
            end
            A(nodes,nodes) = A(nodes,nodes) + ALdiv+ALreg;
    
            % b
            %biotdisplacement
            %modes=[2*nodes(1)-1,2*nodes(1),2*nodes(2)-1,2*nodes(2),2*nodes(3)-1,2*nodes(3)];
            
            fGauss = [f(GaussNodes(1,1),GaussNodes(2,1)),f(GaussNodes(1,2),GaussNodes(2,2)),f(GaussNodes(1,3),GaussNodes(2,3))];
            VectCoeff=[1,vertices(1,1),vertices(1,2);1,vertices(2,1),vertices(2,2 );1,vertices(3,1),vertices(3,2)]\(fVector(nodes));
            fVectGauss = [VectCoeff'*[1;GaussNodes(:,1)],VectCoeff'*[1;GaussNodes(:,2)],VectCoeff'*[1;GaussNodes(:,3)]];
            uprevCoeff=[1,vertices(1,1),vertices(1,2);1,vertices(2,1),vertices(2,2 );1,vertices(3,1),vertices(3,2)]\(u_prev(nodes));
            uVectGauss = [uprevCoeff'*[1;GaussNodes(:,1)],uprevCoeff'*[1;GaussNodes(:,2)],uprevCoeff'*[1;GaussNodes(:,3)]];
           
            BL=zeros(3,1);
            for j= 1:3
                BL(j) = 1/6*det(ARefTri)*sum((fGauss+fVectGauss+M.*uVectGauss).*GaussValues(j,:));
            end
            b(nodes) = b(nodes) + BL;
    
    
    
            %Boundary
            for i=1:length(Dirichlet(:,1))
                A(Dirichlet(i,1),:)=0; A(Dirichlet(i,1),Dirichlet(i,1))=1;
                b(Dirichlet(i,1))=Dirichlet(i,2);
            end
    
            for j = 1:size(Neumann,1)
                nodes = Neumann(j,:);
                vertices = Coordinates(nodes,:);
                E = vertices(2,:) - vertices(1,:);
                b(nodes) = b(nodes) + norm(E)/2 * g(sum(vertices)/2);
            end
        end
        %[L,U]=lu(A);
    
    
    
    u=A\b;
    end
    

    %% Newton
    if scheme == 3
        for k = 1:number_of_elements
            nodes = Elements(k,:);
            vertices = Coordinates(nodes,:);    
     
            %Transformation to refTriangle
            [ARefTri,bRefTri]=RefTriangleMapInv(vertices(:,1),vertices(:,2));
            ARefTriInv = ARefTri^(-1);
            GaussNodes= [ARefTri*[1/6;1/6]+bRefTri,ARefTri*[2/3;1/6]+bRefTri,ARefTri*[1/6;2/3]+bRefTri];
            MVectCoeff=[1,vertices(1,1),vertices(1,2);1,vertices(2,1),vertices(2,2 );1,vertices(3,1),vertices(3,2)]\(bprime(u_prev(nodes)));
            M = [MVectCoeff'*[1;GaussNodes(:,1)],MVectCoeff'*[1;GaussNodes(:,2)],MVectCoeff'*[1;GaussNodes(:,3)]];
            KPVectCoeff=[1,vertices(1,1),vertices(1,2);1,vertices(2,1),vertices(2,2 );1,vertices(3,1),vertices(3,2)]\(kappaprime(u_prev(nodes)));
            KP = [KPVectCoeff'*[1;GaussNodes(:,1)],KPVectCoeff'*[1;GaussNodes(:,2)],KPVectCoeff'*[1;GaussNodes(:,3)]];
            uprevCoeff=[1,vertices(1,1),vertices(1,2);1,vertices(2,1),vertices(2,2 );1,vertices(3,1),vertices(3,2)]\(u_prev(nodes));
            uVectGauss = [uprevCoeff'*[1;GaussNodes(:,1)],uprevCoeff'*[1;GaussNodes(:,2)],uprevCoeff'*[1;GaussNodes(:,3)]];
           
            
            KappaVectCoeff=[1,vertices(1,1),vertices(1,2);1,vertices(2,1),vertices(2,2 );1,vertices(3,1),vertices(3,2)]\(kappa(u_prev(nodes)));
            Kappa = [KappaVectCoeff'*[1;GaussNodes(:,1)],KappaVectCoeff'*[1;GaussNodes(:,2)],KappaVectCoeff'*[1;GaussNodes(:,3)]];

            ALdiv=zeros(3,3);
            for i =1:3
                for j = 1:3
                    ALdiv(i,j) = 1/6*tau*[Coefficients(i,2),Coefficients(i,3)]*(ARefTriInv*ARefTriInv')*[Coefficients(j,2);Coefficients(j,3)]*det(ARefTri)*sum(Kappa);
                end
            end
            
            ALNewt=zeros(3,3);
            for i =1:3
                for j=1:3
                    ALNewt(i,j)=1/6*det(ARefTri)*tau*(uprevCoeff(2)*Coefficients(j,2)+uprevCoeff(3)*Coefficients(j,2))*sum(KP.*GaussValues(i,:));
                end
            end
    
            ALreg=zeros(3,3);
            for i =1:3
                for j=1:3
                    ALreg(i,j)=1/6*det(ARefTri)*sum(M.*GaussValues(i,:).*GaussValues(j,:));
                end
            end
            A(nodes,nodes) = A(nodes,nodes) + ALdiv+ALreg+ALNewt;
    
            % b
            %biotdisplacement
            %modes=[2*nodes(1)-1,2*nodes(1),2*nodes(2)-1,2*nodes(2),2*nodes(3)-1,2*nodes(3)];
            
            fGauss = [f(GaussNodes(1,1),GaussNodes(2,1)),f(GaussNodes(1,2),GaussNodes(2,2)),f(GaussNodes(1,3),GaussNodes(2,3))];
            VectCoeff=[1,vertices(1,1),vertices(1,2);1,vertices(2,1),vertices(2,2 );1,vertices(3,1),vertices(3,2)]\(fVector(nodes));
            fVectGauss = [VectCoeff'*[1;GaussNodes(:,1)],VectCoeff'*[1;GaussNodes(:,2)],VectCoeff'*[1;GaussNodes(:,3)]];
            
            BL=zeros(3,1);
            Bnewt=zeros(3,1);
            for j= 1:3
                BL(j) = 1/6*det(ARefTri)*sum((fGauss+fVectGauss+M.*uVectGauss).*GaussValues(j,:));
                Bnewt(j)=1/6*det(ARefTri)*tau*(uprevCoeff(2)*Coefficients(j,2)+uprevCoeff(3)*Coefficients(j,2))*sum(KP.*uVectGauss);
            end
            b(nodes) = b(nodes) + BL;
    
    
    
            %Boundary
            for i=1:length(Dirichlet(:,1))
                A(Dirichlet(i,1),:)=0; A(Dirichlet(i,1),Dirichlet(i,1))=1;
                b(Dirichlet(i,1))=Dirichlet(i,2);
            end
    
            for j = 1:size(Neumann,1)
                nodes = Neumann(j,:);
                vertices = Coordinates(nodes,:);
                E = vertices(2,:) - vertices(1,:);
                b(nodes) = b(nodes) + norm(E)/2 * g(sum(vertices)/2);
            end
        end
        %[L,U]=lu(A);
    
    

    
    u=A\b;
        
    end

    %% Local L-scheme
    if scheme == 4
        
        for k = 1:number_of_elements
            nodes = Elements(k,:);
            vertices = Coordinates(nodes,:);
            Lb=max(bprime(u_prev(nodes)));
            bm=min(bprime(u_prev(nodes)));
            kappamin=min(kappa(u_prev(nodes)));
            Lk=max(kappaprime(u_prev(nodes)));
            
            
            [ARefTri,bRefTri]=RefTriangleMapInv(vertices(:,1),vertices(:,2));
            ARefTriInv = ARefTri^(-1);
            GaussNodes= [ARefTri*[1/6;1/6]+bRefTri,ARefTri*[2/3;1/6]+bRefTri,ARefTri*[1/6;2/3]+bRefTri];
            KappaVectCoeff=[1,vertices(1,1),vertices(1,2);1,vertices(2,1),vertices(2,2 );1,vertices(3,1),vertices(3,2)]\(kappa(u_prev(nodes)));
            Kappa = [KappaVectCoeff'*[1;GaussNodes(:,1)],KappaVectCoeff'*[1;GaussNodes(:,2)],KappaVectCoeff'*[1;GaussNodes(:,3)]];

            fGauss = [f(GaussNodes(1,1),GaussNodes(2,1)),f(GaussNodes(1,2),GaussNodes(2,2)),f(GaussNodes(1,3),GaussNodes(2,3))];
            VectCoeff=[1,vertices(1,1),vertices(1,2);1,vertices(2,1),vertices(2,2 );1,vertices(3,1),vertices(3,2)]\(fVector(nodes));
            uprevCoeff=[1,vertices(1,1),vertices(1,2);1,vertices(2,1),vertices(2,2 );1,vertices(3,1),vertices(3,2)]\(u_prev(nodes));
            fVectGauss = [VectCoeff'*[1;GaussNodes(:,1)],VectCoeff'*[1;GaussNodes(:,2)],VectCoeff'*[1;GaussNodes(:,3)]];
            uVectGauss = [uprevCoeff'*[1;GaussNodes(:,1)],uprevCoeff'*[1;GaussNodes(:,2)],uprevCoeff'*[1;GaussNodes(:,3)]];
            eta=max(uprevCoeff(2),uprevCoeff(3));
            alpha=Lb*kappamin; beta = 4*kappamin*bm; sigma =2*bm*eta^2*Lk^2*Lb*tau; delta = 2*kappamin^2*tau; xi = tau^2*eta^2*Lk^2*Lb;
            
            gamma=max(0,(sqrt(beta*(alpha^2*beta+alpha*(-beta*delta+2*beta*xi-sigma*delta+delta^2)+(delta-xi)*(sigma*delta-beta*xi)))-alpha*beta+beta*delta-beta*xi)/(beta*delta));
            L=Lb/(2*(1-gamma)-tau*Lb);
            
            
     
            
            ALdiv=zeros(3,3);
            for i =1:3
                for j = 1:3
                    ALdiv(i,j) = 1/6*tau*[Coefficients(i,2),Coefficients(i,3)]*(ARefTriInv*ARefTriInv')*[Coefficients(j,2);Coefficients(j,3)]*det(ARefTri)*sum(Kappa);
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
            BL=zeros(3,1);
           
            for j= 1:3
                BL(j) = 1/6*det(ARefTri)*sum((fGauss+fVectGauss+L*uVectGauss).*GaussValues(j,:));
            end
            b(nodes) = b(nodes) + BL;
    
    
    
            %Boundary
            for i=1:length(Dirichlet(:,1))
                A(Dirichlet(i,1),:)=0; A(Dirichlet(i,1),Dirichlet(i,1))=1;
                b(Dirichlet(i,1))=Dirichlet(i,2);
            end
    
            for j = 1:size(Neumann,1)
                nodes = Neumann(j,:);
                vertices = Coordinates(nodes,:);
                E = vertices(2,:) - vertices(1,:);
                b(nodes) = b(nodes) + norm(E)/2 * g(sum(vertices)/2);
            end
        end
        %[L,U]=lu(A);
    u=A\b;
    end
    end
end
 