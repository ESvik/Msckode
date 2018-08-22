%% Plot Solutions
usq=u(1:2:2*NNnew).*u(2:2:2*NNnew);
unorm=sqrt(usq);
subplot(1,2,1)
trisurf(Elements,X,Y,unorm)