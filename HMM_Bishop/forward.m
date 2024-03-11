% Zachary Loschinskey
% Anna Devor Rotation
% February 2024
% Neurovascular States Project

% Function to calculate the forward probability of an IOHMM with dual alpha
% function defined relationships between the input and output.

% This function uses scaling factors

function [F_norm, Cn] = forward(HMM, X)
    % Function that calculates the normalized forward probabilities for an IOHMM
    % and the correesponding scaling factors

    % Inputs:
    %           HMM struct w/ A, Pi, E
    %           x = 2xN seq of observations

    % Initialize variables
    [n, del] = size(X);  % Length of chain I SWEAR
    k = length(HMM.A(1,:));  % Number of states

    % Starting points of iterations
    F_norm = zeros(n, k); 
    F_current = zeros(1, k);
    Cn = zeros(1, n);

    for state = 1:k
        % Calculate first forward prob
        F_current(state) = HMM.Pi(k) * mvnpdf(X(1,:), HMM.U{state}, HMM.Sigma{state});
    end

    % Normalize and save normalizing factor Cn
    Cn(1) = sum(F_current);
    F_norm(1, :) = F_current ./ Cn(1);

    for t = 2:n  % For each time point 2->n
        for state = 1:k  % For each state 1->k
            F_current = zeros(1, k);
            for prev_state = 1:k
                F_current(prev_state) = F_norm(t-1, prev_state) * HMM.A(prev_state, state);
            end
            F_norm(t, state) = mvnpdf(X(t,:), HMM.U{state}, HMM.Sigma{state}) * sum(F_current);
        end

        % Calculate normalizing factor and normalize
        Cn(t) = sum(F_norm(t, :));
        F_norm(t, :) = F_norm(t, :) ./ Cn(t);
    end
end
