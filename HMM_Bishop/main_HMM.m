% Zachary Loschinskey
% Dr. Brian Depasquale Rotation
% March 2024
% Implementation fo Gaussian HMM on Old Faithful dataset to serve as a
% springboard for finishing NV States GLM HMM

clear;
close all;

%% Import the Old Faithful Dataset of Eruption Lenght vs. Waiting Time
data=readtable("faithful.csv");
X = normalize([data.waiting, data.eruptions]);

%% Initialize the HMM
% A is kxk transition matrix, Pi is 1xk initial prob matrix, E is kx2
% emission matrix of [mean var] for each state
% Initialize the HMM struct
HMM = struct("A", [], "Pi", [], "U", [], "Sigma", []);

% Assign values to A, Pi, U, and Sigma
HMM.A = [0.2, 0.8; 0.5, 0.5];
HMM.Pi = [0.5, 0.5];

% Define cov and mean for each state
U1 = [-1 -1];
U2 = [1 1];
Sigma1 = [1, 0; 0, 1];
Sigma2 = [0.5, 0; 0, 0.5];

% Store covariance matrices in HMM struct
HMM.U = {U1, U2};
HMM.Sigma = {Sigma1, Sigma2};

% Storeage for EM loop
num_iters = 10;
log_lik_stor = zeros(1,num_iters);

%% EM Loop
for chim = 1:num_iters
    % E-Step
    [xi, gamma, log_likelihood]= E_step(HMM,X);
    
    % M-Step
    [HMM] = M_step(HMM, xi, gamma, X);

    % Storage
    log_lik_stor(chim) = log_likelihood;
end

%% Evaluation
figure()
plot(log_lik_stor)
xlabel("EM Iteration")
ylabel("Log Likelihood")
title("Log Likelihood vs. EM Iteration")

figure()
plot(X(:,1), X(:,2), 'k.')
hold on;
plot_gaus(HMM.U{1}, HMM.Sigma{1})
plot_gaus(HMM.U{2}, HMM.Sigma{2})
xlabel("Time Between Eruptions (mins)")
ylabel("Eruption Time (mins)")
title("Old Faithful Eruption Time vs. Waiting Time")

