p1=plot(analysis(3:end-1,1),analysis(3:end-1,2),'b*-')
hold on
p2=plot(analysis(2:end-1,3),analysis(2:end-1,4),'ro-')
p3=plot(analysis(2:end-1,5),analysis(2:end-1,6),'k:')
p4=plot(analysis(2:end-1,7),analysis(2:end-1,8),'g--')
p5=plot(analysis(end,1),analysis(end,2),'p','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10)
p6=plot(analysis(end,3),analysis(end,4),'p','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10)
p7=plot(analysis(end,5),analysis(end,6),'p','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',10)
p8=plot(analysis(end,7),analysis(end,8),'p','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',10)
hold off
legend([p1,p2,p3,p4],'\kappa=0.01','\kappa=0.1','\kappa=1','\kappa=10')
ylabel('Iterations')
xlabel('L/L_b')
set(gca,'fontsize',15)