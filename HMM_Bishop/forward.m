% Zachary Loschinskey
% Anna Devor Rotation
% February 2024
% Neurovascular States Project

% Function to calculate the forward probability of an IOHMM with dual alpha
% function defined relationships between the input and output.

% This function uses scaling factors

function [alpha, c] = forward(HMM, X)
    % Function that calculates the normalized forward probabilities for an IOHMM
    % and the correesponding scaling factors

    % Inputs:
    %           HMM struct w/ A, Pi, E
    %           x = 2xN seq of observations

    % Initialize variables
    [n, ~] = size(X);  % Length of chain I SWEAR
    k = length(HMM.A(1,:));  % Number of states

    % Starting points of iterations
    alpha = zeros(n, k); 
    c = zeros(1, n);

    for state = 1:k
        % Calculate first forward prob
        alpha(1, state) = HMM.Pi(state) * mvnpdf(X(1,:), HMM.U{state}, HMM.Sigma{state});
    end

    % Normalize and save normalizing factor Cn
    c(1) = 1 / (sum(alpha(1, :)) + eps);
    alpha(1, :) = alpha(1, :) * c(1);

    for t = 2:n  % For each time point 2->n
        for state = 1:k  % For each state 1->k
            alpha(t, state) = 0.0;
            for prev_state = 1:k
                alpha(t, state) = alpha(t, state) + alpha(t-1, prev_state) * HMM.A(prev_state, state);
            end
            alpha(t, state) = alpha(t, state) * mvnpdf(X(t,:), HMM.U{state}, HMM.Sigma{state});
        end

        % Calculate normalizing factor and normalize
        c(t) = 1 / (sum(alpha(t, :)) + eps);
       alpha(t, :) = alpha(t, :) * c(t);
    end
end
