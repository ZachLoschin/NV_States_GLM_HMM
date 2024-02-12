% Zachary Loschinskey
% Anna Devor Rotation
% January 2024
% Neurovascular States Project

% Function to calculate the E_Step for fitting alpha function data
% Synthetic data has been created to test this function in cooperation
% with the M_step_alpha()

% When finished, these functions should be able to be implemented for real
% NV state learning on Ca and Hb data windowed at 7s sampled at 10Hz.


function [xi_11, xi_12, xi_21, xi_22, gamma1, gamma2] = E_step_alpha(input, output, A, E, Pi, latent)
    
    % Define the number of states
    num_states = length(A(:,1));    
    
    % Define how many timepoints there are, each is 70 datapoints long
    T = length(input(:, 1));
    disp("Calculating Alpha")

    % Calculate the forwad and backward probabilities of the IOHMM
    [alpha_norm, Cn] = forward_iterative(input, output, A, E, Pi);  % Checking if i need to switch input and output

    disp("Calculating Beta")
    beta_norm = backward_iterative(input, output, A, E, Cn);

    % -- Interestingly the beta_norm values barely don't row-wise sum to 1
    
    % Calculate the gamma values for each class
    gamma1 = alpha_norm(:,1) .* beta_norm(:,1);
    gamma2 = alpha_norm(:,2) .* beta_norm(:,2);

    n = length(input);

    % Initialize xi matrices
    xi_11 = zeros(n-1, 1);
    xi_12 = zeros(n-1, 1);
    xi_21 = zeros(n-1, 1);
    xi_22 = zeros(n-1, 1);
    
    % Calculate xi values using the formula in Bishop
    for t = 1:n-1
        xi_11(t) = alpha_norm(t, 1) * A(1, 1) * alpha_emit_prob(input(t+1, :), output(t+1, :), E(1,:)) * beta_norm(t+1, 1) / Cn(t);
        xi_12(t) = alpha_norm(t, 1) * A(1, 2) * alpha_emit_prob(input(t+1, :), output(t+1, :), E(1,:)) * beta_norm(t+1, 2) / Cn(t);
        xi_21(t) = alpha_norm(t, 2) * A(2, 1) * alpha_emit_prob(input(t+1, :), output(t+1, :), E(2,:)) * beta_norm(t+1, 1) / Cn(t);
        xi_22(t) = alpha_norm(t, 2) * A(2, 2) * alpha_emit_prob(input(t+1, :), output(t+1, :), E(2,:)) * beta_norm(t+1, 2) / Cn(t);
    end

        

end
