% find the average fitness cost for all four motifs
% 1st entry: non-autoregulated positive control
% 2nd entry: non-autoregulated negative control
% 3rd entry:     autoregulated positive control
% 4th entry:     autoregulated negative control

% in case of regulatory delay, x needs to be shifted, e.g., using circshift

function fitness = find_fitness(nu_m,nu_p,s_p,s_r,D,T)
        
        % discretize time
        tau = T/1e4;        
        ti = [tau:tau:T];
              
        % positive control without self-activation                   
        T_s = D*T;
        T_n = (1-D)*T;
        s = s_p; % for solving x
        s_b  = s_r*(ti<=T_n) + s_r*(ti>T_n);
        s_nb = s_r*(ti<=T_n) + (s_r+s_p)*(ti>T_n);
        x = find_x(ti,nu_m,nu_p,s,T_n,T_s);
        fitness(1) = 1/T*trapz(ti,x.*s_nb + (1-x).*s_b); 


        % negative control without self-repression
        T_s = (1-D)*T;
        T_n = D*T;
        s = s_p; % for solving x
        s_b  = s_r*(ti<=T_n) + s_r*(ti>T_n);
        s_nb = s_r*(ti<=T_n) + (s_r+s_p)*(ti>T_n);        
        x = find_x(ti,nu_m,nu_p,s,T_n,T_s);
        fitness(2) = 1/T*trapz(ti,x.*s_nb + (1-x).*s_b); 
         
        
        % positive control with self-activation
        T_s = D*T;
        T_n = (1-D)*T;
        s = s_p - s_r; % for solving x        
        s_b  = 0*(ti<=T_n) + s_r*(ti>T_n);
        s_nb = 0*(ti<=T_n) + (0+s_p)*(ti>T_n);
        x = find_x(ti,nu_m,nu_p,s,T_n,T_s);
        fitness(3) = 1/T*trapz(ti,x.*s_nb + (1-x).*s_b); 
        

        % negative control with self-repression
        T_s = (1-D)*T;
        T_n = D*T;
        s = s_p + s_r; % for solving x
        s_b  = s_r*(ti<=T_n) + 0*(ti>T_n);
        s_nb = s_r*(ti<=T_n) + (s_r+s_p)*(ti>T_n);
        x = find_x(ti,nu_m,nu_p,s,T_n,T_s);    
        fitness(4) = 1/T*trapz(ti,x.*s_nb + (1-x).*s_b); 
end