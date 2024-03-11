% Zachary Loschinskey
% Dr. Anna Devor Rotation
% January 2024
% Training of IOHMM and Evaluation on the test dataset
% EM Implementation of fitting alpha function datasets

% Uses log of values rather than scaling factors

clear;
clc;
close all;

%% Dataset Generation

% Set the seed for reproducibility
rng(1);

% Generate the alpha dataset
[Ca, Hb, latent, trans, E_real] = gen_linear_data();



%% Initialization of Model Parameters and Storage Cell Arrays
% Initialize the parameters A, E, and Pi
% Transition matrix initialization --> Playing with this shows weird LL
% behavior...
A = [0.7 .3;
    0.2 0.8];

% Initial guess dual alpha function parameters 
% m1 b1; m2 b2
E = [5 1;
    8 0.5];

% Generate initial guess dual alpha functions 1x70 for each class
x = 1:1:100;

% Define the two linear relationships for initial guess
y1 = E(1,1).*x + E(1,2);
y2 = E(2,1).*x + E(2,2);

% Initial state distribution initialization
Pi = [0.5 0.5];

% Storage variables for viewing training history
PiCell = {};
ACell = {};
ECell = {};
weight_diff = [];  % -- Average difference in weights btwn class1 and class2
log_likelihood_storage = [];
state_prediction_storage = [];
mean_g1 = [];
mean_g2 = [];
weights1 = [];
weights2 = [];

[N, del] = size(Ca);

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
    [xi, gamma, alpha_log, beta_log] = E_step_linear(Ca, Hb, A, E, Pi, latent);

    log_likelihood = logsumexp(alpha_log(N, :));
    %--------------------------------------------------------------------------
    % M-Step: Maximize the Q(theta_old, theta) function
    %         w/ respect to pi, A, E, and U
    %--------------------------------------------------------------------------
    [Pi_new, A_new, E_new, weightDiff, state_prediction, w1, w2] = M_step_linear(Ca, Hb, xi, gamma, E);

    %--------------------------------------------------------------------------
    % Save history of parameters and update variables for next iteration
    %--------------------------------------------------------------------------
    Pi = Pi_new;
    A = A_new;
    E = E_new;
    weight_diff(chim) = weightDiff;
    log_likelihood_storage = [log_likelihood_storage, log_likelihood];
    state_prediction_storage = [state_prediction_storage, state_prediction];
    %mean_g1 = [mean_g1, mean(gamma1)];
    %mean_g2 = [mean_g2, mean(gamma2)];
    weights1 = [weights1, w1];
    weights2 = [weights2, w2];

end

% -- Normalize the log-likelihood and diff values -- %
log_likelihood_storage = log_likelihood_storage(2:end);
norm_LL = rescale(log_likelihood_storage);
d_LL = diff(log_likelihood_storage);
norm_diff = rescale(d_LL);

% -- Plot the fits over time -- %
plot_fit_linear(ECell, E_real)

% -- Plot the log likelihood over time -- %
figure()
plot(log_likelihood_storage)
title("Log Likelihood Vs. EM Iterations")
ylabel("Log Likelihood")
xlabel("EM Iterations")
legend("LL")

% -- Final Predictions of States -- %
final_pred = state_prediction_storage(:,end);

state_prediction_storage = state_prediction_storage == 2;

% -- Plot the transition matrix predicted values over time -- %
% Initialize vectors to store extracted elements
A_11 = [];
A_12 = [];
A_21 = [];
A_22 = [];

% Loop through each cell in ACell
for i = 1:numel(ACell)
    % Extract the 2x2 double from the current cell
    current_A = ACell{i};
    
    % Extract elements and store them
    A_11 = [A_11, current_A(1,1)];
    A_12 = [A_12, current_A(1,2)];
    A_21 = [A_21, current_A(2,1)];
    A_22 = [A_22, current_A(2,2)];
end

% -- Calculate the real transition matrix from the data -- %
prev = 1;
a_data = [0 0;
          0 0];

for idx = 2:length(latent)
    current = latent(idx);

    if current == 1 && prev == 1
        a_data(1, 1) = a_data(1,1) + 1;
    end

    if current == 1 && prev == 2
        a_data(2, 1) = a_data(2,1) + 1;
    end

    if current == 2 && prev == 1
        a_data(1, 2) = a_data(1,2) + 1;
    end

    if current == 2 && prev == 2
        a_data(2,2) = a_data(2,2) + 1;
    end

    prev = current;

end


a_data = a_data ./ sum(a_data,2);

A_real_11 = repelem(a_data(1,1), length(ACell));
A_real_22 = repelem(a_data(2,2), length(ACell));

% Plot all elements on the same plot with different colors
figure;
hold on;
plot(A_11, 'r', 'DisplayName', 'A(1,1)');
%plot(A_12, 'r--', 'DisplayName', 'A(1,2)');
%plot(A_21, 'b--', 'DisplayName', 'A(2,1)');
plot(A_22, 'b', 'DisplayName', 'A(2,2)');
plot(A_real_11, 'r*' ,"DisplayName", 'Real A(1,1)');
plot(A_real_22, 'b*',"DisplayName", 'Real A(2,2)');
title("Transition probabilities over EM iteration")
hold off;
ylim([0 1])

% Add legend
legend('show', "Location", "best");










