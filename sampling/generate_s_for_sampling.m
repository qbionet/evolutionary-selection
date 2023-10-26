%% generate s(t) for all variants
% first row contains s for the non-binders
% second row contains s for the binders
function s_for_sampling = generate_s_for_sampling(variant, T_ni, T_i, s_p, s_r, a, cycles)
    switch variant

        % activation, no autoregulation
        case '+no'
            T_n = T_ni;
            T_s = T_i;
            s_b  = [s_r*ones(1,T_n) s_r*ones(1,T_s)];
            s_nb = [s_r*ones(1,T_n) (s_r+s_p)*ones(1,T_s)];

        % repression, no autoregulation    
        case '-no'
            T_n = T_i;
            T_s = T_ni;
            s_b  = [s_r*ones(1,T_n) s_r*ones(1,T_s)];
            s_nb = [s_r*ones(1,T_n) (s_r+s_p)*ones(1,T_s)];

        % activation, with autoregulation    
        case '+fb'
            T_n = T_ni;
            T_s = T_i;
            s_b  = [0*ones(1,T_n) a*s_r*ones(1,T_s)];
            s_nb = [0*ones(1,T_n) (0+s_p)*ones(1,T_s)];

        % repression, with autoregulation    
        case '-fb'
            T_n = T_i;
            T_s = T_ni;
            s_b  = [a*s_r*ones(1,T_n) 0*ones(1,T_s)];
            s_nb = [a*s_r*ones(1,T_n) (a*s_r+s_p)*ones(1,T_s)];
    end

    s(1,:) = repmat(s_nb,1,cycles); % non-binder
    s(2,:) = repmat(s_b,1,cycles);  % binder
    
    s_for_sampling = s(1,:)-s(2,:);
    
end