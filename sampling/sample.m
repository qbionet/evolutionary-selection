% stochastic sampling using a Wright-Fisher model
function n = Sample(N,nu_m,nu_p,s,n0)
    sat = @(s) min(max(s, 0), 1);
    
    n(:,1) = n0;
    for t = 1:length(s)
        n_current = n(:,t);
        s_current = s(t);
    
        % step 1
        N_nb = n_current(1); % non-binder
        N_b  = n_current(2); % binder

        m_m = poissrnd(N_b*nu_m);   % b->nb mutation
        m_p = poissrnd(N_nb*nu_p);  % bn->b mutation

        x = sat((N_nb + m_m - m_p)/N);
            
        % step 2
        x_new = sat(x - s_current*x*(1-x)/(1 - s_current*x));
        n(1,t+1) = binornd(N,x_new);
        n(2,t+1) = N - n(1,t+1);
    end
    n(:,1) = [];

end