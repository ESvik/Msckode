%% Edge numbering Mandel

function [Elements,Coordinates] = ElementsPlusEdgesMandel(Elements,Coordinates,dx,dy,xmin,xmax,ymin,ymax)
NN=((xmax-xmin)/dx+1)*((ymax-ymin)/dy+1);
Nx=(xmax-xmin)/dx+1; Ny=(ymax-ymin)/dy+1;
i=Elements(:,1);
j=Elements(:,2);
k=Elements(:,3);
A=sparse(i,j,-1,NN,NN);
A=A+sparse(i,k,-1,NN,NN);
A=A+sparse(j,k,-1,NN,NN);
A=A+A';
A=triu(A);
[r,c,v]=find(A);
Entries=1:length(v);
Entries=Entries+NN;
A=sparse(r,c,Entries,NN,NN);
A=A+A';
NewElements=zeros(length(Elements(:,1)),6);
for i =1:length(Elements(:,1))
    NewElements(i,:)=[Elements(i,1),Elements(i,2),Elements(i,3),A(Elements(i,1),Elements(i,2)),A(Elements(i,1),Elements(i,3)),A(Elements(i,2),Elements(i,3))];
end
Elements=NewElements;
NewCoordinates=zeros((2*Nx-1)*(2*Ny-1),2);
NewCoordinates(1:length(Coordinates),:)=Coordinates;
for i = 1:length(Elements(:,1))
    Nodes=Elements(i,:);
    NewCoordinates(Nodes(4),:)=0.5*(Coordinates(Nodes(1),:)+Coordinates(Nodes(2),:));
    NewCoordinates(Nodes(5),:)=0.5*(Coordinates(Nodes(1),:)+Coordinates(Nodes(3),:));
    NewCoordinates(Nodes(6),:)=0.5*(Coordinates(Nodes(2),:)+Coordinates(Nodes(3),:));
end
Coordinates=NewCoordinates;
end