% Code created by Loïc Marrec

% Parameters

n = 2e0;                        % Hill coefficient
theta = 1e3;                    % Inflection time
gW = 0.1;                       % Death rate of W microbes
XW_i = 10;                      % Initial number of W microbes
fG = .5;                        % Reproduction rate of G microbes 
gG = 0.1;                       % Death rate of G microbes
XG_i = 0;                       % Initial number of G microbes
K = 1e3;                        % Carrying capacity
t0 = 50;                        % Mutant appearance time
Nit = 1e3;                      % Number of stochastic realizations  

% Simulation

pfix = Gillespie_fct(Nit, n, theta, gW, XW_i, fG, gG, XG_i, K, t0);
