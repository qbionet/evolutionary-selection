% Compare average fitness cost in the deterministic case
clear all, close all, clc

nu_m = 1e-7;      % loss-of-function mutation rate
nu_p = nu_m/10;   % gain-of-function mutation rate
s_p  = 10*nu_m;   % maximum P-cost
D    = 0.95;      % demand

T_vector = logspace(0,10,10);
control_vector = linspace(0,1,10);
[T_grid, control_grid] = meshgrid(T_vector,control_vector);
for i = 1:size(T_grid,1)
    for j = 1:size(T_grid,2)
        T   = T_grid(i,j);              % period length
        s_r =s_p*control_grid(i,j);     % maximum R-cost
        
        % average fitness cost over a period
        % 1st entry: non-autoregulated positive control
        % 2nd entry: non-autoregulated negative control
        % 3rd entry:     autoregulated positive control
        % 4th entry:     autoregulated negative control
        fitness = find_fitness(nu_m,nu_p,s_p,s_r,D,T);
        
        F(i,j,:) = fitness;
    end
end
