% Diffusion approximation
clear all, close all, clc


nu_m = 1e-3;    % rate of loss-of-function mutations
nu_p = nu_m/10; % rate of gain-of-function mutations

D = 0.05;        % demand
T = 1e4;         % period length
s_p = 100*nu_m;  % max P-cost
s_r = 0.1*s_p;   % max R-cost
N = 1000;        % population size
Ne = N/2;        % effective population size


% find the average fitness cost for all four motifs
% 1st entry: non-autoregulated positive control
% 2nd entry: non-autoregulated negative control
% 3rd entry:     autoregulated positive control
% 4th entry:     autoregulated negative control
fitness = find_fitness(nu_m,nu_p,s_p,s_r,D,T,Ne)