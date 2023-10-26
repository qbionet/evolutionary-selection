% average fitness cost for all four motifs using stochastic sampling
clear all, close all, clc

T   = 100;   % period length
N   = 1e3;   % population size

nu_m = 1e-3;     % rate of loss-of-function mutations
nu_p = 0.1*nu_m; % rate of gain-of-function mutations
s_p = 100*nu_m;  % max P-cost
s_r = 0.25*s_p;  % max R-cost
a   = 2;         % factor increase in the R-cost due to autoregulation

sampling_rep = 100;    % number of sampling repetitions
sampling_cycles = 100; % number of sampling cycles

% low demand case
D   = 0.05;
base_variant = '+'; 
[x_pos_noFB x_pos_FB fitness_pos_noFB fitness_pos_FB] = find_fitness(nu_m, nu_p, T, N, D, s_p, s_r, a, sampling_rep, sampling_cycles, base_variant);

% high demand case
D   = 0.95;
base_variant = '-'; 
[x_neg_noFB x_neg_FB fitness_neg_noFB fitness_neg_FB] = find_fitness(nu_m, nu_p, T, N, D, s_p, s_r, a, sampling_rep, sampling_cycles, base_variant);