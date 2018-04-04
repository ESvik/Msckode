%Theoretical analysis

function rate= theoreticalraterichards(gamma,L_b,b_m,C,k,tau)

rate = (L_b-4*(1-gamma).*b_m.*gamma)./(L_b+4*(1-gamma)*(tau*k)/C);

end