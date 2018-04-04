%gives the transformation from the triangle (xi,yi), (xj,yj), (xk,yk) to
%the triangle (0,0),(1,0),(0,1)

function [A,b]=RefTriangleMap(x,y)
    E=[x(1),y(1),1;x(2),y(2),1;x(3),y(3),1];
    A=zeros(2,2);
    b=zeros(2,1);
    u1=E\[0;1;0];
    u2=E\[0;0;1];
    A(1,:)=[u1(1),u1(2)]; A(2,:)=[u2(1),u2(2)]; b=[u1(3);u2(3)];
end