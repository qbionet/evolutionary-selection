% constraint for boundary conditions at 0/T and Tn
function F = AB_constraint_fun(z,s,nu_m,nu,s_tilde,s_hat,Ts,Tn)
    F(1) = nu_m/nu + z(1)*exp(-nu*Tn) - 1/(2*s)*(s_tilde - s_hat*(1-z(2))/(1+z(2)));
    F(2) = nu_m/nu + z(1) - 1/(2*s)*(s_tilde - s_hat*(1-z(2)*exp(-s_hat*Ts))/(1+z(2)*exp(-s_hat*Ts)));
end