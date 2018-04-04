
%Maps the ref triangle to x(1),y(1) x(2),y(2) x(3),y(3)

function [A,b]=RefTriangleMapInv(x,y)
    b=[x(1);y(1)];
    A=[x(2)-x(1),x(3)-x(1);y(2)-y(1),y(3)-y(1)];
end