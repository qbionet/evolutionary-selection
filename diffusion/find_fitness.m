% find the average fitness cost for all four motifs
% 1st entry: non-autoregulated positive control
% 2nd entry: non-autoregulated negative control
% 3rd entry:     autoregulated positive control
% 4th entry:     autoregulated negative control
function fitness = find_fitness(nu_m,nu_p,s_p,s_r,D,T,Ne)

    N_samplingT = 1e3;

    tau = T/N_samplingT;
    K = 30; % spatial grid for solving the PDE numerically

    ti = tau:tau:T;
    nu   = nu_m + nu_p;

    x0 = 0.5;   % initial starting point
    f0 = [zeros(1,round(x0*K)) 1 zeros(1,K-round(x0*K))]'*K; % initial probability

    xi = linspace(0,1,K+1);
    max_cycle = 10+1;

    iterations = 1:max_cycle*N_samplingT;
    last_cycle_index = (max_cycle-1)*N_samplingT+1:length(iterations);
    time       = iterations*tau;
    max_iter   = length(iterations);



    % + control without feedback
    T_s = D*T;
    T_n = T-T_s;
    s_b  = s_r*(ti<=T_n) + s_r*(ti>T_n);
    s_nb = s_r*(ti<=T_n) + (s_r+s_p)*(ti>T_n);

    s_for_x = s_p;
    s_cycle   = [repmat(0,1,round(T_n/tau)) repmat(s_for_x,1,round(T_s/tau))];
    s_t       = repmat(s_cycle,1,max_cycle);

    S   = repmat(s_t,K+1,1);
    X   = repmat(xi',1,max_iter);
    M = nu_m - nu*X - S.*X.*(1-X);
    [f P] = Waxman(f0,K,tau,Ne,max_iter,M,last_cycle_index);
    x = find_x(ti,nu_m,nu_p,s_for_x,T_n,T_s);
    for k = 1:length(x)
        x_bar(k) = P(end,k) + trapz(xi,xi'.*P(:,k));
    end                
    fitness(1) = 1/T*trapz(ti,x_bar.*s_nb + (1-x_bar).*s_b); 
    


    % - control without feedback
    T_n = D*T;
    T_s = T-T_n;
    s_b  = s_r*(ti<=T_n) + s_r*(ti>T_n);
    s_nb = s_r*(ti<=T_n) + (s_r+s_p)*(ti>T_n);

    s_for_x = s_p;
    s_cycle   = [repmat(0,1,round(T_n/tau)) repmat(s_for_x,1,round(T_s/tau))];
    s_t       = repmat(s_cycle,1,max_cycle);
    
    S   = repmat(s_t,K+1,1);
    X   = repmat(xi',1,max_iter);
    M = nu_m - nu*X - S.*X.*(1-X);
    [f P] = Waxman(f0,K,tau,Ne,max_iter,M,last_cycle_index);
    x = find_x(ti,nu_m,nu_p,s_for_x,T_n,T_s);
    for k = 1:length(x)
        x_bar(k) = P(end,k) + trapz(xi,xi'.*P(:,k));
    end                
    fitness(2)  = 1/T*trapz(ti,x_bar.*s_nb + (1-x_bar).*s_b); 

    

    % + control with feedback
    T_s = D*T;
    T_n = T-T_s;
    s_b  = 0*(ti<=T_n) + s_r*(ti>T_n);
    s_nb = 0*(ti<=T_n) + (0+s_p)*(ti>T_n);

    s_for_x = s_p-s_r;
    s_cycle   = [repmat(0,1,round(T_n/tau)) repmat(s_for_x,1,round(T_s/tau))];
    s_t       = repmat(s_cycle,1,max_cycle);
    
    S   = repmat(s_t,K+1,1);
    X   = repmat(xi',1,max_iter);
    M = nu_m - nu*X - S.*X.*(1-X);
    [f P] = Waxman(f0,K,tau,Ne,max_iter,M,last_cycle_index);
    x = find_x(ti,nu_m,nu_p,s_for_x,T_n,T_s);
    for k = 1:length(x)
        x_bar(k) = P(end,k) + trapz(xi,xi'.*P(:,k));
    end                
    fitness(3)  = 1/T*trapz(ti,x_bar.*s_nb + (1-x_bar).*s_b); 
    


    % - control with feedback
    T_n = D*T;
    T_s = T-T_n;
    s_b  = s_r*(ti<=T_n) + 0*(ti>T_n);
    s_nb = s_r*(ti<=T_n) + (s_r+s_p)*(ti>T_n);

    s_for_x = s_p+s_r;
    s_cycle   = [repmat(0,1,round(T_n/tau)) repmat(s_for_x,1,round(T_s/tau))];
    s_t       = repmat(s_cycle,1,max_cycle);
    
    S   = repmat(s_t,K+1,1);
    X   = repmat(xi',1,max_iter);
    M = nu_m - nu*X - S.*X.*(1-X);
    [f P] = Waxman(f0,K,tau,Ne,max_iter,M,last_cycle_index);
    x = find_x(ti,nu_m,nu_p,s_for_x,T_n,T_s);
    for k = 1:length(x)
        x_bar(k) = P(end,k) + trapz(xi,xi'.*P(:,k));
    end                
    fitness(4)  = 1/T*trapz(ti,x_bar.*s_nb + (1-x_bar).*s_b); 
end