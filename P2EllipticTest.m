%% P2 Elliptic test
f = @(x,y) -2*x.*(x-1)-2*y.*(y-1);
uexact = @(x,y) (x-1).*x.*y.*(y-1);
x_min=0; x_max=1; y_min=0; y_max=1;
Error = [];
for h=[1/2,1/4,1/8,1/16,1/32]
[x,y]=meshgrid(x_min:h:x_max,y_min:h:y_max);
X=reshape(x,[],1); Y=reshape(y,[],1);
DT=delaunayTriangulation(X,Y);
Coordinates=DT.Points;
Elements=DT.ConnectivityList;
Neumann=[];
g=0;
NN=length(X);
[ElementsP2,CoordinatesP2]=ElementsPlusEdges(Elements,Coordinates,NN);

%Basis and its evaluation points for P2
CoefficientsP2= zeros(6,6);
MatOfCoord=[1,0,0,0,0,0;1,1,0,1,0,0;1,0,1,0,1,0;1,1/2,0,1/4,0,0;1,0,1/2,0,1/4,0;1,1/2,1/2,1/4,1/4,1/4];   
CoefficientsP2(1,:)=MatOfCoord\[1;0;0;0;0;0];
CoefficientsP2(2,:)=MatOfCoord\[0;1;0;0;0;0];
CoefficientsP2(3,:)=MatOfCoord\[0;0;1;0;0;0];
CoefficientsP2(4,:)=MatOfCoord\[0;0;0;1;0;0];
CoefficientsP2(5,:)=MatOfCoord\[0;0;0;0;1;0];
CoefficientsP2(6,:)=MatOfCoord\[0;0;0;0;0;1];
GaussValuesP2 = zeros(6,3);
for i = 1:6
    GaussValuesP2(i,1) = (CoefficientsP2(i,1) + CoefficientsP2(i,2)*1/6 + CoefficientsP2(i,3)*1/6+CoefficientsP2(i,4)*1/6^2+CoefficientsP2(i,5)*1/6^2+CoefficientsP2(i,6)*1/6^2);
    GaussValuesP2(i,2) = (CoefficientsP2(i,1) + CoefficientsP2(i,2)*2/3 + CoefficientsP2(i,3)*1/6+CoefficientsP2(i,4)*(2/3)^2+CoefficientsP2(i,5)*1/6^2+CoefficientsP2(i,6)*1/6*2/3);
    GaussValuesP2(i,3) = (CoefficientsP2(i,1) + CoefficientsP2(i,2)*1/6 + CoefficientsP2(i,3)*2/3+CoefficientsP2(i,4)*1/6^2+CoefficientsP2(i,5)*(2/3)^2+CoefficientsP2(i,6)*1/6*2/3);
end
Dirichletu = zeros(4*(2/h+1),2);
Dirichletu(:,1)=[find(CoordinatesP2(:,1)==x_min);find(CoordinatesP2(:,1)==x_max);find(CoordinatesP2(:,2)==y_min);find(CoordinatesP2(:,2)==y_max)];
Dirichletu(1:(2/h+1),2)=0;
Dirichletu((2/h+1)+1:2*(2/h+1),2)=0;
Dirichletu(2*(2/h+1)+1:3*(2/h+1),2)=0;
Dirichletu(3*(2/h+1)+1:4*(2/h+1),2)=0;



u=FEMElliptic2DP2(CoordinatesP2,ElementsP2,1,Dirichletu,f,CoefficientsP2,GaussValuesP2,NN);

Error=[Error,norm(u-uexact(CoordinatesP2(:,1),CoordinatesP2(:,2)),2)]
end

%trisurf(Elements,CoordinatesP2(:,1),CoordinatesP2(:,2),u)
