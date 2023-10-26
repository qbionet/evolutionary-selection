%% find x
function x = find_x(t,nu_m,nu_p,s,Tn,Ts)
    nu = nu_p + nu_m;
    if s
        s_tilde = s + nu;
        s_hat = sqrt(s_tilde^2- 4*s*nu_m);
        
        fun = @(z)AB_constraint_fun(z,s,nu_m,nu,s_tilde,s_hat,Ts,Tn);
        options = optimset('Display','off');
        z = fsolve(fun,[0 0],options);
        A = z(1);
        B = z(2);
            
        xn = nu_m/nu + A*exp(-nu*t);
        xs = 1/(2*s)*(s_tilde - s_hat*(1-B*exp(-s_hat*(t-Tn)))./(1+B*exp(-s_hat*(t-Tn))));
        
        x = [xn(find(t<=Tn)) xs(find(t>Tn))];
    else
        x = repmat(nu_m/(nu_m+nu_p),1,length(t));
    end
end