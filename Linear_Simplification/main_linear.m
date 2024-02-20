% Zachary Loschinskey
% Dr. Brian Depasquale
% January 2024
% Training of IOHMM and Evaluation on the test dataset
% EM Implementation of fitting linear functions

clear;
clc;
close all;

%% Dataset Generation

% Set the seed for reproducibility
rng(1);

% Generate the alpha dataset
[Ca, Hb, latent, y1, y2, A_real, E_real] = gen_linear_data();

% Plot the alpha functions
E_gen = E_real;


%% Initialization of Model Parameters and Storage Cell Arrays
% Initialize the parameters A, E, and Pi
% Transition matrix initialization --> Playing with this shows weird LL
% behavior...
A = [0.6 0.4;
    0.29 0.71];

% Initial guess dual alpha function parameters 
% m1 b1 ; m2 b2
E = [5 0;
    -2 0.5];

% x vector for fitting
x = 0.1:0.1:7;

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


%% Expectation-Maximization Algorithm
for chim = 1:100
    PiCell{chim} = Pi;
    ACell{chim} = A;
    ECell{chim} = E;

    % disp(chim)
    %--------------------------------------------------------------------------
    % E-Step: Calculate P(Z|X,Theta)
    %         Includes calculating gamma(Zn) and eta(Zn-1, Zn)
    %--------------------------------------------------------------------------
    [xi_11, xi_12, xi_21, xi_22, gamma1, gamma2, log_likelihood] = E_step_linear(Ca, Hb, A, E, Pi, latent);
    
    %--------------------------------------------------------------------------
    % M-Step: Maximize the Q(theta_old, theta) function
    %         w/ respect to pi, A, E, and U
    %--------------------------------------------------------------------------
    [Pi_new, A_new, E_new, gamma1, gamma2, w1, w2] = M_step_linear(Ca, Hb, xi_11, xi_12, xi_21, xi_22, gamma1, gamma2, E);

    %--------------------------------------------------------------------------
    % Save history of parameters and update variables for next iteration
    %--------------------------------------------------------------------------
    Pi = Pi_new;
    A = A_new;
    E = E_new;
    log_likelihood_storage = [log_likelihood_storage, log_likelihood];
    mean_g1 = [mean_g1, mean(gamma1)];
    mean_g2 = [mean_g2, mean(gamma2)];
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

% -- Plot normalized LL and Diff values -- %
figure()
plot(norm_LL)
hold on;
plot(norm_diff)
title("Normalized LL and Change in LL")
xlabel("EM Iterations")
legend(["Normalized LL", "Normalized Diff LL"], "Location","best")
ylim([-0.05 1.05])


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

A_real_11 = A_real(1,1);
A_real_22 = A_real(2,2);

A_real_11 = repelem(A_real_11, length(ACell));
A_real_22 = repelem(A_real_22, length(ACell));

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



%% Analysis
% % Initialize a matrix to store the differences
% differences1 = zeros(size(weights1, 1), size(weights1, 2) - 1);
% 
% % Compute differences between subsequent columns
% for i = 1:size(weights1, 2) - 1
%     differences1(:, i) = weights1(:, i+1) - weights1(:, i);
% end
% 
% 
% % Initialize a matrix to store the differences
% differences2 = zeros(size(weights2, 1), size(weights2, 2) - 1);
% 
% % Compute differences between subsequent columns
% for i = 1:size(weights2, 2) - 1
%     differences2(:, i) = weights2(:, i+1) - weights2(:, i);
% end
% 
% d1 = mean(differences1);
% d2 = mean(differences2);
% figure()
% hold on
% plot(d1)
% plot(d2)
% title("Mean Change in IRF for Each Class")
% ylabel("Mean Change in Weights")
% xlabel("EM Iteration")
% legend("Class1", "Class2")




