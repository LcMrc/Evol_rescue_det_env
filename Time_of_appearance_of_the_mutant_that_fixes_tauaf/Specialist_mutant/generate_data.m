% Code created by Loïc Marrec

% Parameters

n = 5;                          % Hill coefficient
theta = 1e3;                    % Inflection time
gW = 0.1;                       % Death rate of W microbes
XW_i = 10;                      % Initial number of W microbes
gS = 0.1;                       % Death rate of G microbes
XS_i = 0;                       % Initial number of G microbes
K = 1e3;                        % Carrying capacity
mu = 1e-5;                      % Mutation probability upon division
Nit = 1e3;                      % Number of stochastic realizations 

% Simulation

tauaf = Gillespie_fct(Nit, n, theta, gW, XW_i, gS, XS_i, K, mu);
