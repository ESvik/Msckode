%% Getting A, D and M_p 

function [A,D,M_p]=BiotSchur(Coordinates,Elements,ElementsP2,lambda,mu,Dirichletp,Dirichletu)
NN=length(Coordinates(:,1));
number_of_elements=length(Elements(:,1));

%Basis and its evaluation points for P1
Coefficients= zeros(3,3);
MatOfCoord=[1,0,0;1,1,0;1,0,1];
for i =1:3
    bb=zeros(3,1);
    bb(i)=1;
    Coefficients(i,:)=MatOfCoord\bb;
end
GaussValues = zeros(3,3);
for i = 1:3
    GaussValues(i,1) = 1/6*(Coefficients(i,1) + Coefficients(i,2)*1/6 + Coefficients(i,3)*1/6);
    GaussValues(i,2) = 1/6*(Coefficients(i,1) + Coefficients(i,2)*2/3 + Coefficients(i,3)*1/6);
    GaussValues(i,3) = 1/6*(Coefficients(i,1) + Coefficients(i,2)*1/6 + Coefficients(i,3)*2/3);
end
    %Basis and its evaluation points for P2
CoefficientsP2= zeros(6,6);
MatOfCoord=[1,0,0,0,0,0;1,1,0,1,0,0;1,0,1,0,1,0;1,1/2,0,1/4,0,0;1,0,1/2,0,1/4,0;1,1/2,1/2,1/4,1/4,1/4];   
CoefficientsP2(1,:)=MatOfCoord\[1;0;0;0;0;0];
CoefficientsP2(2,:)=MatOfCoord\[0;1;0;0;0;0];
CoefficientsP2(3,:)=MatOfCoord\[0;0;1;0;0;0];
CoefficientsP2(4,:)=MatOfCoord\[0;0;0;1;0;0];
CoefficientsP2(5,:)=MatOfCoord\[0;0;0;0;1;0];
CoefficientsP2(6,:)=MatOfCoord\[0;0;0;0;0;1];
CoefficientsP2=[CoefficientsP2(1,:);CoefficientsP2(1,:);CoefficientsP2(2,:);CoefficientsP2(2,:);CoefficientsP2(3,:);CoefficientsP2(3,:);CoefficientsP2(4,:);CoefficientsP2(4,:);CoefficientsP2(5,:);CoefficientsP2(5,:);CoefficientsP2(6,:);CoefficientsP2(6,:)];
GaussValuesP2 = zeros(12,3);
for i = 1:12
    GaussValuesP2(i,1) = (CoefficientsP2(i,1) + CoefficientsP2(i,2)*1/6 + CoefficientsP2(i,3)*1/6+CoefficientsP2(i,4)*1/6^2+CoefficientsP2(i,5)*1/6^2+CoefficientsP2(i,6)*1/6^2);
    GaussValuesP2(i,2) = (CoefficientsP2(i,1) + CoefficientsP2(i,2)*2/3 + CoefficientsP2(i,3)*1/6+CoefficientsP2(i,4)*(2/3)^2+CoefficientsP2(i,5)*1/6^2+CoefficientsP2(i,6)*1/6*2/3);
    GaussValuesP2(i,3) = (CoefficientsP2(i,1) + CoefficientsP2(i,2)*1/6 + CoefficientsP2(i,3)*2/3+CoefficientsP2(i,4)*1/6^2+CoefficientsP2(i,5)*(2/3)^2+CoefficientsP2(i,6)*1/6*2/3);
end

%Assembly
A=sparse(2*(2*sqrt(NN)-1)^2,2*(2*sqrt(NN)-1)^2);
D = sparse(2*(2*sqrt(NN)-1)^2,NN);
M_p = sparse(NN,NN);

for k = 1:number_of_elements
        for k = 1:number_of_elements
            nodes = ElementsP2(k,:);
            originalnodes=[nodes(1),nodes(2),nodes(3)];
            vertices = Coordinates(originalnodes,:);
            nodesp2=[2*nodes(1)-1,2*nodes(1),2*nodes(2)-1,2*nodes(2),2*nodes(3)-1,2*nodes(3),2*nodes(4)-1,2*nodes(4),2*nodes(5)-1,2*nodes(5),2*nodes(6)-1,2*nodes(6)];
            nodes = Elements(k,:);
               
            %Transformation to refTriangle
            [ARefTri,~]=RefTriangleMapInv(vertices(:,1),vertices(:,2));
            [ARefTriinv,~]=RefTriangleMap(vertices(:,1),vertices(:,2));
            %ARefTriinv=ARefTri^(-1);
            
            
            %A
            dphix = @(x,y,i) ARefTriinv(1,1)*(CoefficientsP2(i,2)+2*CoefficientsP2(i,4)*x+CoefficientsP2(i,6)*y)+ARefTriinv(2,1)*(CoefficientsP2(i,3)+2*CoefficientsP2(i,5)*y+CoefficientsP2(i,6)*x);
            dphiy = @(x,y,i) ARefTriinv(1,2)*(CoefficientsP2(i,2)+2*CoefficientsP2(i,4)*x+CoefficientsP2(i,6)*y)+ARefTriinv(2,2)*(CoefficientsP2(i,3)+2*CoefficientsP2(i,5)*y+CoefficientsP2(i,6)*x);

            
            ALreg=zeros(12,12);
            for i =1:2:11
                for j = 1:2:11
                    ip = @(x,y) dphix(x,y,i)*dphix(x,y,j) + 1/2*dphiy(x,y,i)*dphiy(x,y,j);
                    ALreg(i,j) = 2/6*mu*det(ARefTri)*(ip(1/6,1/6)+ip(1/6,2/3)+ip(2/3,1/6));
                    %ALdiv(i,j)=1/6*B*det(ARefTri)*sum([(ARefTriinv(1,1)*(Coefficients(i,2)+2*Coefficients(i,4)*1/6+Coefficients(i,6)*1/6)+ARefTriinv(2,1)*(Coefficients(i,3)+2*Coefficients(i,5)*1/6+Coefficients(i,6)*1/6))*(ARefTriinv(1,1)*(Coefficients(j,2)+2*Coefficients(j,4)*1/6+Coefficients(j,6)*1/6)+ARefTriinv(2,1)*(Coefficients(j,3)+2*Coefficients(j,5)*1/6+Coefficients(j,6)*1/6)),(ARefTriinv(1,1)*(Coefficients(i,2)+2*Coefficients(i,4)*1/6+Coefficients(i,6)*2/3)+ARefTriinv(2,1)*(Coefficients(i,3)+2*Coefficients(i,5)*2/3+Coefficients(i,6)*1/6))*(ARefTriinv(1,1)*(Coefficients(j,2)+2*Coefficients(j,4)*1/6+Coefficients(j,6)*2/3)+ARefTriinv(2,1)*(Coefficients(j,3)+2*Coefficients(j,5)*2/3+Coefficients(j,6)*1/6)),(ARefTriinv(1,1)*(Coefficients(i,2)+2*Coefficients(i,4)*2/3+Coefficients(i,6)*1/6)+ARefTriinv(2,1)*(Coefficients(i,3)+2*Coefficients(i,5)*1/6+Coefficients(i,6)*2/3))*(ARefTriinv(1,1)*(Coefficients(j,2)+2*Coefficients(j,4)*2/3+Coefficients(j,6)*1/6)+ARefTriinv(2,1)*(Coefficients(j,3)+2*Coefficients(j,5)*1/6+Coefficients(j,6)*2/3))]);
                end
            end
            for i =1:2:11
                for j = 2:2:12
                    ip = @(x,y) 1/2*dphiy(x,y,i)*dphix(x,y,j);
                    ALreg(i,j) = 2/6*mu*det(ARefTri)*(ip(1/6,1/6)+ip(1/6,2/3)+ip(2/3,1/6));
                    %ALdiv(i,j)=1/6*B*det(ARefTri)*sum([(ARefTriinv(1,1)*(Coefficients(i,2)+2*Coefficients(i,4)*1/6+Coefficients(i,6)*1/6)+ARefTriinv(2,1)*(Coefficients(i,3)+2*Coefficients(i,5)*1/6+Coefficients(i,6)*1/6))*(ARefTriinv(1,2)*(Coefficients(j,2)+2*Coefficients(j,4)*1/6+Coefficients(j,6)*1/6)+ARefTriinv(2,2)*(Coefficients(j,3)+2*Coefficients(j,5)*1/6+Coefficients(j,6)*1/6)),(ARefTriinv(1,1)*(Coefficients(i,2)+2*Coefficients(i,4)*1/6+Coefficients(i,6)*2/3)+ARefTriinv(2,1)*(Coefficients(i,3)+2*Coefficients(i,5)*2/3+Coefficients(i,6)*1/6))*(ARefTriinv(1,2)*(Coefficients(j,2)+2*Coefficients(j,4)*1/6+Coefficients(j,6)*2/3)+ARefTriinv(2,2)*(Coefficients(j,3)+2*Coefficients(j,5)*2/3+Coefficients(j,6)*1/6)),(ARefTriinv(1,1)*(Coefficients(i,2)+2*Coefficients(i,4)*2/3+Coefficients(i,6)*1/6)+ARefTriinv(2,1)*(Coefficients(i,3)+2*Coefficients(i,5)*1/6+Coefficients(i,6)*2/3))*(ARefTriinv(1,2)*(Coefficients(j,2)+2*Coefficients(j,4)*2/3+Coefficients(j,6)*1/6)+ARefTriinv(2,2)*(Coefficients(j,3)+2*Coefficients(j,5)*1/6+Coefficients(j,6)*2/3))]);    
                end
            end
            for i =2:2:12
                for j = 1:2:11
                    ip = @(x,y) 1/2*dphix(x,y,i)*dphiy(x,y,j);
                    ALreg(i,j) = 2/6*mu*det(ARefTri)*(ip(1/6,1/6)+ip(1/6,2/3)+ip(2/3,1/6));
                    %ALdiv(i,j)=1/6*B*det(ARefTri)*sum([(ARefTriinv(1,2)*(Coefficients(i,2)+2*Coefficients(i,4)*1/6+Coefficients(i,6)*1/6)+ARefTriinv(2,2)*(Coefficients(i,3)+2*Coefficients(i,5)*1/6+Coefficients(i,6)*1/6))*(ARefTriinv(1,1)*(Coefficients(j,2)+2*Coefficients(j,4)*1/6+Coefficients(j,6)*1/6)+ARefTriinv(2,1)*(Coefficients(j,3)+2*Coefficients(j,5)*1/6+Coefficients(j,6)*1/6)),(ARefTriinv(1,2)*(Coefficients(i,2)+2*Coefficients(i,4)*1/6+Coefficients(i,6)*2/3)+ARefTriinv(2,2)*(Coefficients(i,3)+2*Coefficients(i,5)*2/3+Coefficients(i,6)*1/6))*(ARefTriinv(1,1)*(Coefficients(j,2)+2*Coefficients(j,4)*1/6+Coefficients(j,6)*2/3)+ARefTriinv(2,1)*(Coefficients(j,3)+2*Coefficients(j,5)*2/3+Coefficients(j,6)*1/6)),(ARefTriinv(1,2)*(Coefficients(i,2)+2*Coefficients(i,4)*2/3+Coefficients(i,6)*1/6)+ARefTriinv(2,2)*(Coefficients(i,3)+2*Coefficients(i,5)*1/6+Coefficients(i,6)*2/3))*(ARefTriinv(1,1)*(Coefficients(j,2)+2*Coefficients(j,4)*2/3+Coefficients(j,6)*1/6)+ARefTriinv(2,1)*(Coefficients(j,3)+2*Coefficients(j,5)*1/6+Coefficients(j,6)*2/3))]);
                end
            end
            for i =2:2:12
                for j = 2:2:12
                    ip = @(x,y) dphiy(x,y,i)*dphiy(x,y,j) + 1/2*dphix(x,y,i)*dphix(x,y,j);
                    ALreg(i,j) = 2/6*mu*det(ARefTri)*(ip(1/6,1/6)+ip(1/6,2/3)+ip(2/3,1/6));
                    %ALdiv(i,j)=1/6*B*det(ARefTri)*sum([(ARefTriinv(1,2)*(Coefficients(i,2)+2*Coefficients(i,4)*1/6+Coefficients(i,6)*1/6)+ARefTriinv(2,2)*(Coefficients(i,3)+2*Coefficients(i,5)*1/6+Coefficients(i,6)*1/6))*(ARefTriinv(1,2)*(Coefficients(j,2)+2*Coefficients(j,4)*1/6+Coefficients(j,6)*1/6)+ARefTriinv(2,2)*(Coefficients(j,3)+2*Coefficients(j,5)*1/6+Coefficients(j,6)*1/6)),(ARefTriinv(1,2)*(Coefficients(i,2)+2*Coefficients(i,4)*1/6+Coefficients(i,6)*2/3)+ARefTriinv(2,2)*(Coefficients(i,3)+2*Coefficients(i,5)*2/3+Coefficients(i,6)*1/6))*(ARefTriinv(1,2)*(Coefficients(j,2)+2*Coefficients(j,4)*1/6+Coefficients(j,6)*2/3)+ARefTriinv(2,2)*(Coefficients(j,3)+2*Coefficients(j,5)*2/3+Coefficients(j,6)*1/6)),(ARefTriinv(1,2)*(Coefficients(i,2)+2*Coefficients(i,4)*2/3+Coefficients(i,6)*1/6)+ARefTriinv(2,2)*(Coefficients(i,3)+2*Coefficients(i,5)*1/6+Coefficients(i,6)*2/3))*(ARefTriinv(1,2)*(Coefficients(j,2)+2*Coefficients(j,4)*2/3+Coefficients(j,6)*1/6)+ARefTriinv(2,2)*(Coefficients(j,3)+2*Coefficients(j,5)*1/6+Coefficients(j,6)*2/3))]);
                end
            end
     
            ALdiv=zeros(12,12);
            for i =1:2:11
                for j=1:2:11
                    ip = @(x,y) dphix(x,y,i)*dphix(x,y,j);
                    ALdiv(i,j) = 1/6*lambda*det(ARefTri)*(ip(1/6,1/6)+ip(1/6,2/3)+ip(2/3,1/6));
                    %ALdiv(i,j)=1/6*det(ARefTri)*sum([(ARefTriinv(1,1)*dphix(1/6,1/6,i)+ARefTriinv(2,1)*dphiy(1/6,1/6,i))*(ARefTriinv(1,1)*dphix(1/6,1/6,j)+ARefTriinv(2,1)*dphiy(1/6,1/6,j))+1/2*(dphiy(1/6,1/6,i)*ARefTriinv(2,2)+dphix(1/6,1/6,i)*ARefTriinv(1,2))*(dphiy(1/6,1/6,j)*ARefTriinv(2,2)+dphix(1/6,1/6,j)*ARefTriinv(1,2)),(ARefTriinv(1,1)*dphix(2/3,1/6,i)+ARefTriinv(2,1)*dphiy(2/3,1/6,i))*(ARefTriinv(1,1)*dphix(2/3,1/6,j)+ARefTriinv(2,1)*dphiy(2/3,1/6,j))+1/2*(dphiy(2/3,1/6,i)*ARefTriinv(2,2)+dphix(2/3,1/6,i)*ARefTriinv(1,2))*(dphiy(2/3,1/6,j)*ARefTriinv(2,2)+dphix(2/3,1/6,j)*ARefTriinv(1,2)),(ARefTriinv(1,1)*dphix(1/6,2/3,i)+ARefTriinv(2,1)*dphiy(1/6,2/3,i))*(ARefTriinv(1,1)*dphix(1/6,2/3,j)+ARefTriinv(2,1)*dphiy(1/6,2/3,j))+1/2*(dphiy(1/6,2/3,i)*ARefTriinv(2,2)+dphix(1/6,2/3,i)*ARefTriinv(1,2))*(dphiy(1/6,2/3,j)*ARefTriinv(2,2)+dphix(1/6,2/3,j)*ARefTriinv(1,2))]);
                end
            end
            for i =2:2:12
                for j=2:2:12
                    ip = @(x,y) dphiy(x,y,i)*dphiy(x,y,j);
                    ALdiv(i,j) = 1/6*lambda*det(ARefTri)*(ip(1/6,1/6)+ip(1/6,2/3)+ip(2/3,1/6));
                    %ALdiv(i,j)=1/6*det(ARefTri)*sum([(dphix(1/6,1/6,i)*ARefTriinv(1,2)+dphiy(1/6,1/6,i)*ARefTriinv(2,2))*(dphix(1/6,1/6,j)*ARefTriinv(1,2)+dphiy(1/6,1/6,j)*ARefTriinv(2,2))+1/2*(dphiy(1/6,1/6,i)*ARefTriinv(2,1)+dphix(1/6,1/6,i)*ARefTriinv(1,1))*(dphiy(1/6,1/6,j)*ARefTriinv(2,1)+dphix(1/6,1/6,j)*ARefTriinv(1,1)),(dphix(2/3,1/6,i)*ARefTriinv(1,2)+dphiy(2/3,1/6,i)*ARefTriinv(2,2))*(dphix(2/3,1/6,j)*ARefTriinv(1,2)+dphiy(2/3,1/6,j)*ARefTriinv(2,2))+1/2*(dphiy(2/3,1/6,i)*ARefTriinv(2,1)+dphix(2/3,1/6,i)*ARefTriinv(1,1))*(dphiy(2/3,1/6,j)*ARefTriinv(2,1)+dphix(2/3,1/6,j)*ARefTriinv(1,1)),(dphix(1/6,2/3,i)*ARefTriinv(1,2)+dphiy(1/6,2/3,i)*ARefTriinv(2,2))*(dphix(1/6,2/3,j)*ARefTriinv(1,2)+dphiy(1/6,2/3,j)*ARefTriinv(2,2))+1/2*(dphiy(1/6,2/3,i)*ARefTriinv(2,1)+dphix(1/6,2/3,i)*ARefTriinv(1,1))*(dphiy(1/6,2/3,j)*ARefTriinv(2,1)+dphix(1/6,2/3,j)*ARefTriinv(1,1))]);
                end
            end
            for i =1:2:11
                for j=2:2:12
                    ip = @(x,y) dphix(x,y,i)*dphiy(x,y,j);
                    ALdiv(i,j) = 1/6*lambda*det(ARefTri)*(ip(1/6,1/6)+ip(1/6,2/3)+ip(2/3,1/6));
                    %ALdiv(i,j)=1/12*det(ARefTri)*sum([(dphiy(1/6,1/6,i)*ARefTriinv(2,2)+dphix(1/6,1/6,i)*ARefTriinv(1,2))*(dphiy(1/6,1/6,j)*ARefTriinv(2,1)+dphix(1/6,1/6,j)*ARefTriinv(1,1)),(dphiy(2/3,1/6,i)*ARefTriinv(2,2)+dphix(2/3,1/6,i)*ARefTriinv(1,2))*(dphiy(2/3,1/6,j)*ARefTriinv(2,1)+dphix(2/3,1/6,j)*ARefTriinv(1,1)),(dphiy(1/6,2/3,i)*ARefTriinv(2,2)+dphix(1/6,2/3,i)*ARefTriinv(1,2))*(dphiy(1/6,2/3,j)*ARefTriinv(2,1)+dphix(1/6,2/3,j)*ARefTriinv(1,1))]);
                end
            end
            for i =2:2:12
                for j=1:2:11
                    ip = @(x,y) dphiy(x,y,i)*dphix(x,y,j);
                    ALdiv(i,j) = 1/6*lambda*det(ARefTri)*(ip(1/6,1/6)+ip(1/6,2/3)+ip(2/3,1/6));
                    %ALdiv(i,j)=1/12*det(ARefTri)*sum([(dphiy(1/6,1/6,j)*ARefTriinv(2,2)+dphix(1/6,1/6,j)*ARefTriinv(1,2))*(dphiy(1/6,1/6,i)*ARefTriinv(2,1)+dphix(1/6,1/6,i)*ARefTriinv(1,1)),(dphiy(2/3,1/6,j)*ARefTriinv(2,2)+dphix(2/3,1/6,j)*ARefTriinv(1,2))*(dphiy(2/3,1/6,i)*ARefTriinv(2,1)+dphix(2/3,1/6,i)*ARefTriinv(1,1)),(dphiy(1/6,2/3,j)*ARefTriinv(2,2)+dphix(1/6,2/3,j)*ARefTriinv(1,2))*(dphiy(1/6,2/3,i)*ARefTriinv(2,1)+dphix(1/6,2/3,i)*ARefTriinv(1,1))]);
                end
            end

            A(nodesp2,nodesp2) = A(nodesp2,nodesp2) + ALdiv+ALreg;
            
            %M_p
            ALreg=zeros(3,3);
            for i =1:3
                for j=1:3
                    ALreg(i,j)=6*det(ARefTri)*sum(GaussValues(i,:).*GaussValues(j,:));
                end
            end
            M_p(nodes,nodes) = M_p(nodes,nodes) +ALreg;
            
            %D
            Dreg=zeros(12,3);
            for i =1:3
                for j=1:2:11
                    Dreg(j,i) = det(ARefTri)*sum(GaussValues(i,:).*[(ARefTriinv(1,1)*(CoefficientsP2(i,2)+2*CoefficientsP2(i,4)*1/6+CoefficientsP2(i,6)*1/6)+ARefTriinv(2,1)*(CoefficientsP2(i,3)+2*CoefficientsP2(i,5)*1/6+CoefficientsP2(i,6)*1/6)),(ARefTriinv(1,1)*(CoefficientsP2(i,2)+2*CoefficientsP2(i,4)*2/3+CoefficientsP2(i,6)*1/6)+ARefTriinv(2,1)*(CoefficientsP2(i,3)+2*CoefficientsP2(i,5)*2/3+CoefficientsP2(i,6)*1/6)),(ARefTriinv(1,1)*(CoefficientsP2(i,2)+2*CoefficientsP2(i,4)*1/6+CoefficientsP2(i,6)*2/3)+ARefTriinv(2,1)*(CoefficientsP2(i,3)+2*CoefficientsP2(i,5)*1/6+CoefficientsP2(i,6)*2/3))]);
                end
            end
            for i =1:3
                for j=2:2:12
                    Dreg(j,i) = det(ARefTri)*sum(GaussValues(i,:).*[(ARefTriinv(1,2)*(CoefficientsP2(i,2)+2*CoefficientsP2(i,4)*1/6+CoefficientsP2(i,6)*1/6)+ARefTriinv(2,2)*(CoefficientsP2(i,3)+2*CoefficientsP2(i,5)*1/6+CoefficientsP2(i,6)*1/6)),(ARefTriinv(1,2)*(CoefficientsP2(i,2)+2*CoefficientsP2(i,4)*2/3+CoefficientsP2(i,6)*1/6)+ARefTriinv(2,2)*(CoefficientsP2(i,3)+2*CoefficientsP2(i,5)*2/3+CoefficientsP2(i,6)*1/6)),(ARefTriinv(1,2)*(CoefficientsP2(i,2)+2*CoefficientsP2(i,4)*1/6+CoefficientsP2(i,6)*2/3)+ARefTriinv(2,2)*(CoefficientsP2(i,3)+2*CoefficientsP2(i,5)*1/6+CoefficientsP2(i,6)*2/3))]);
                end
            end
            D(nodesp2,nodes)=D(nodesp2,nodes)+Dreg;
end


for i=1:length(Dirichletu)
    A(Dirichletu(i,1),:)=0; A(Dirichletu(i,1),Dirichletu(i,1))=1;
    D(Dirichletu(i,1),:)=0;
end
for i=1:length(Dirichletp)
    M_p(Dirichletp(i),:)=0; M_p(Dirichletp(i),Dirichletp(i))=1;
    D(:,Dirichletp(i))=0;
end


end
           