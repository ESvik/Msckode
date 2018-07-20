%%PLot generator

p1 = plot(Analysis(1:end-1,1),Analysis(1:end-1,2))
hold on
p2 = plot(Analysis(1:end-1,3),Analysis(1:end-1,4),'-*')
p3 = plot(Analysis(1:end-1,5),Analysis(1:end-1,6),'-o')
p4 = plot(Analysis(1:end-1,7),Analysis(1:end-1,8),'.:')
p5 = plot(Analysis(1:end-1,9),Analysis(1:end-1,10),'-.')
p6 = plot(Analysis(1:end-1,11),Analysis(1:end-1,12),'--')
hold off
legend([p1,p2,p3,p4,p5,p6],'\kappa = 10^{-15}','\kappa = 10^{-14}','\kappa = 10^{-13}','\kappa = 10^{-12}','\kappa = 10^{-11}','\kappa = 10^{-10}')
xlabel('\delta')
ylabel('Iterations')
title('Analysis of scheme with different \kappa')