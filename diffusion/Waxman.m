%% Waxman-solver
function [f P] = Waxman(f0,K,tau,Ne,max_iter,M,last_cycle_index)
    
    xi = linspace(0,1,K+1);
    epsilon = 1/K;
    alpha = tau/(2*epsilon^2);

    f(:,1) = f0;    % initial probability distribution

    for iteration = 1:max_iter
        M_actual = M(:,iteration)';
        U = epsilon*M_actual/2 - xi.*(1-xi)/(4*Ne);
        V = epsilon*M_actual/2 + xi.*(1-xi)/(4*Ne);
        
        U_mod = U(2:K+1);
        U_mod(1) = 2*U_mod(1);
        V_mod = V(1:K);
        V_mod(K) = 2*V_mod(K);
        
        R = diag(V-U) + diag(U_mod,1) - diag(V_mod,-1);
        R(1,1) = 2*V(1);
        R(K+1,K+1) = -2*U(K+1);
    
        f_current = f(:,iteration);
        f_next = inv(eye(K+1) + alpha*R)*(eye(K+1) - alpha*R)*f_current;
        f(:,iteration+1) = f_next;
    end

    f(:,1) = [];
    P = f(:,last_cycle_index)*epsilon;
    P(1,:) = P(1,:)/2;
    P(end,:) = P(end,:)/2;

end