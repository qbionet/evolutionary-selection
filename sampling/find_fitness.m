function [x_noFB x_FB s_noFB s_FB] = find_fitness(nu_m, nu_p, T, N, D, s_p, s_r, a, sampling_rep, sampling_cycles, base_variant)
    x0 = 0.5;   % initial starting point
    sampling_ti     = 1:1:T;

    edges = [0:0.01:1]; % histogram bin edges
    bincenters = (edges(1:end-1) + edges(2:end))/2;

    drop_cycles = 3; % drop first 3 cycles
    keep_index  = drop_cycles*T+1:sampling_cycles*T; % keep the rest

    T_i = round(D*T);
    T_ni = T-T_i;

    switch base_variant
        case '-'
            T_s = T_ni;
            T_n = T_i;

            sb_noFB  = s_r*(sampling_ti<=T_n) + s_r*(sampling_ti>T_n);
            snb_noFB = s_r*(sampling_ti<=T_n) + (s_r+s_p)*(sampling_ti>T_n);

            sb_FB  = a*s_r*(sampling_ti<=T_n) + 0*(sampling_ti>T_n);
            snb_FB = a*s_r*(sampling_ti<=T_n) + (a*s_r+s_p)*(sampling_ti>T_n);

        case '+'
            T_s = T_i;
            T_n = T_ni;

            sb_noFB  = s_r*(sampling_ti<=T_n) + s_r*(sampling_ti>T_n);
            snb_noFB = s_r*(sampling_ti<=T_n) + (s_r+s_p)*(sampling_ti>T_n);

            sb_FB  = 0*(sampling_ti<=T_n) + a*s_r*(sampling_ti>T_n);
            snb_FB = 0*(sampling_ti<=T_n) + (0+s_p)*(sampling_ti>T_n);
    end

    % control without feedback
    variant = [base_variant 'no']; 
    s_for_sampling = generate_s_for_sampling(variant, T_ni, T_i, s_p, s_r, a, sampling_cycles);       
    n0 = [round(0.5*N) N-round(0.5*N)];
        
    sb_int  = repmat(sb_noFB,1,sampling_cycles-drop_cycles);
    snb_int = repmat(snb_noFB,1,sampling_cycles-drop_cycles);
        
    for rep_index = 1:sampling_rep
        n = sample(N,nu_m,nu_p,s_for_sampling,n0);
        x_rep(rep_index,:) = n(1,:)/N;
        x_rep_truncated(rep_index,:) = x_rep(rep_index,keep_index);

        x_int   = squeeze(x_rep_truncated(rep_index,:));
        s_rep_noFB(rep_index)  = sum(x_int.*snb_int + (1-x_int).*sb_int)/length(keep_index);        
    end
        
    data = squeeze(mean(x_rep_truncated(:,:)));
    hcount = histcounts(data(:),edges);
    x_noFB(1,:) = hcount;
    s_noFB = mean(s_rep_noFB);

    % control with feedback
    variant = [base_variant 'fb']; 
    s_for_sampling = generate_s_for_sampling(variant, T_ni, T_i, s_p, s_r, a, sampling_cycles);       
    n0 = [round(0.5*N) N-round(0.5*N)];
        
    sb_int  = repmat(sb_FB,1,sampling_cycles-drop_cycles);
    snb_int = repmat(snb_FB,1,sampling_cycles-drop_cycles);
        
    for rep_index = 1:sampling_rep
        n = sample(N,nu_m,nu_p,s_for_sampling,n0);
        x_rep(rep_index,:) = n(1,:)/N;
        x_rep_truncated(rep_index,:) = x_rep(rep_index,keep_index);

        x_int   = squeeze(x_rep_truncated(rep_index,:));        
        s_rep_FB(rep_index)  = sum(x_int.*snb_int + (1-x_int).*sb_int)/length(keep_index);        
    end
        
    data = squeeze(mean(x_rep_truncated(:,:)));
    hcount = histcounts(data(:),edges);
    x_FB(1,:) = hcount;
    s_FB = mean(s_rep_FB);    
end