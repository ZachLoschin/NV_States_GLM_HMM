% Zachary Loschinskey
% Dr. Anna Devor Rotation
% January 2024
% Training of IOHMM and Evaluation on the test dataset
% EM Implementation of fitting alpha function datasets

clear;
clc;
close all;

%% Dataset Generation
% Generate the alpha dataset
[Ca, Hb, latent, alpha1, alpha2] = gen_alpha_data();

% Plot the alpha functions
E_real = [alpha1; alpha2];


%% Initialization of Model Parameters and Storage Cell Arrays
% Initialize the parameters A, E, and Pi
A = [0.4 0.6;
    0.45 0.55];


% Trans probabilities of synthetic data
% transition_probabilities = [0.75, 0.25;
                            % 0.12, 0.88]

% Initial alpha parameters
% E = [0.4 0.5;
%     0.3 0.4];

E = [0.7 0.75;
    0.8 0.9];

% E = [0.5 0.7;
%     0.1 0.8];

x = 0.1:0.1:7;

alpha1_t0 = 0;
alpha1_tau1 = E(1,1);
alpha1_tau2 = E(1,2);

alpha1_1 = (((x-alpha1_t0) ./ alpha1_tau1) .^3) .* exp(-(x-alpha1_t0) ./ alpha1_tau1);
alpha1_2 = -(((x-alpha1_t0) ./ alpha1_tau2) .^3) .* exp(-(x-alpha1_t0) ./ alpha1_tau2);
alpha1 = alpha1_1 + alpha1_2;


alpha1_t0 = 0;
alpha1_tau1 = E(2,1);
alpha1_tau2 = E(2,2);

alpha1_1 = (((x-alpha1_t0) ./ alpha1_tau1) .^3) .* exp(-(x-alpha1_t0) ./ alpha1_tau1);
alpha1_2 = -(((x-alpha1_t0) ./ alpha1_tau2) .^3) .* exp(-(x-alpha1_t0) ./ alpha1_tau2);
alpha2 = alpha1_1 + alpha1_2;


E = [alpha1;alpha2];

Pi = [0.5 0.5];

% Storage variables for training history
PiCell = {};
ACell = {};
ECell = {};
weight_diff = [];  % Store the max weight difference
log_likelihood_storage = [];

%% Expectation-Maximization Algorithm
for chim = 1:50
    PiCell{chim} = Pi;
    ACell{chim} = A;
    ECell{chim} = E;

    disp(chim)
    %--------------------------------------------------------------------------
    % E-Step: Calculate P(Z|X,Theta)
    %         Includes calculating gamma(Zn) and eta(Zn-1, Zn)
    %--------------------------------------------------------------------------
    [xi_11, xi_12, xi_21, xi_22, gamma1, gamma2, log_likelihood] = E_step_alpha(Ca, Hb, A, E, Pi, latent);
    
    %--------------------------------------------------------------------------
    % M-Step: Maximize the Q(theta_old, theta) function 
    %         w/ respect to pi, A, E, and U  
    %--------------------------------------------------------------------------
    [Pi_new, A_new, E_new, weightDiff] = M_step_alpha(Ca, Hb, xi_11, xi_12, xi_21, xi_22, gamma1, gamma2, E);

    %--------------------------------------------------------------------------
    % Save history of parameters and update variables for next iteration
    %--------------------------------------------------------------------------
    
    Pi = Pi_new;
    A = A_new;
    E = E_new;
    weight_diff(chim) = weightDiff;
    log_likelihood_storage = [log_likelihood_storage, log_likelihood];
    
end

% -- Plot the fits over time -- %
plot_fit2(ECell, E_real)

% -- Plot the log likelihood over time -- %
figure()
plot(log_likelihood_storage)
title("Log Likelihood Vs. EM Iterations")
ylabel("Log Likelihood")
xlabel("EM Iterations")











