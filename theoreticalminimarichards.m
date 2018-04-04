%Calculating gammas
function gamma = theoreticalminimarichards(L_b,b_m,C,tau,k)
gamma=1+(L_b*b_m*C-sqrt(L_b^2*b_m^2*C^2+4*L_b*b_m^2*C*k*tau+4*L_b*b_m*k^2*tau^2))/(4*b_m*k*tau);
end