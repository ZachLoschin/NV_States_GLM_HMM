% Zachary Loschinskey
% Anna Devor Rotation
% January 2024
% Neurovascular States Project

% Function to calculate the E_Step for fitting alpha function data
% Synthetic data has been created to test this function in cooperation
% % with the M_step_alpha()

% When finished, these functions should be able to be implemented for real
% NV state learning on Ca and Hb data windowed at 7s sampled at 10Hz.


function [xi, gamma, alpha_log, beta_log] = E_step_linear(input, output, A, E, Pi, latent)
    
    % Define the number of states
    num_states = length(A(:,1));    
    
    % Define how many timepoints there are, each is 70 datapoints long
    T = length(input(:, 1));

    % Calculate the forwad and backward probabilities of the IOHMM
    alpha_log = forward_iterative_linear(input, output, A, E, Pi);

    beta_log = backward_iterative_linear(input, output, A, E);

    gamma = alpha_log + beta_log;  % Each column is for a state
    
    % -- Normalize the gamma values -- %
    for t = 1:length(gamma)
        % Calculate the maximum value of gamma at time step t
        max_gamma = max(gamma(t, :));

        % Calculate the normalization factor (log-sum-exp trick)
        log_sum = max_gamma + log(sum(exp(gamma(t, :) - max_gamma)));

        % Subtract log_sum from each element of gamma(:, t)
        normalized_gamma(t, :) = gamma(t, :) - log_sum;
    end
    gamma = normalized_gamma;
    [N, del] = size(input);
    [K, del] = size(A);
    
    xi = zeros(K,K,N-1);
    

    for t = 1:N-1
        log_xi_unnormalized = zeros(K, K);

        for i = 1:K
            for j = 1:K
                                                                                                                        % WRONG
                log_xi_unnormalized(i, j) = alpha_log(t, i) + log(A(i,j)) + linear_emit(input(t+1, :), output(t+1, :), E(K,:)) + beta_log(t+1,j);
            end
        end

        % Normalize the xi values
        max_xi = max(max(log_xi_unnormalized));
        denom = max_xi + log(sum(sum(exp(log_xi_unnormalized - max_xi))));
        xi(:,:,t) = log_xi_unnormalized - denom;    % Note that values of the bottom row are always a factor of 10 larger than the top row
    end

end
    


