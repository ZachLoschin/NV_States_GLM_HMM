% Zachary Loschinskey
% Anna Devor Rotation
% February 2024
% Neurovascular States Project

% Function to calculate the forward probability of an IOHMM with dual alpha
% function defined relationships between the input and output.

% This function uses scaling factors

function [F_norm, Cn] = forward_iterative(x, u, A, E, Pi)
    % Function that calculates the normalized forward probabilities for an IOHMM
    % and the correesponding scaling factors

    % Inputs:   x = Nx70 vector of observed outputs
    %           u = Nx70 vector of inputs
    %           A = KxK matrix of transition probabilities
    %           E = 2xK matrix of emission probabilities

    % Initialize variables
    n = length(x);  % Length of chain
    k = length(A(1,:));  % Number of states

    % Starting points of iterations
    F_norm = zeros(n, k); 
    F_current = zeros(1, k);
    Cn = zeros(1, n);

    for state = 1:k
        % Calculate first forward prob
        F_current(state) = Pi(k) * alpha_emit_prob(x(1, :), u(1, :), E(state,:));
    end

    % Normalize and save normalizing factor Cn
    Cn(1) = sum(F_current);
    F_norm(1, :) = F_current ./ Cn(1);

    for t = 2:n  % For each time point 2->n
        for state = 1:k  % For each state 1->k
            F_current = zeros(1, k);
            for prev_state = 1:k
                F_current(prev_state) = F_norm(t-1, prev_state) * A(prev_state, state);
            end
            F_norm(t, state) = alpha_emit_prob(x(t, :), u(t, :), E(state, :)) * sum(F_current);
        end

        % Calculate normalizing factor and normalize
        Cn(t) = sum(F_norm(t, :));
        F_norm(t, :) = F_norm(t, :) ./ Cn(t);
    end
end
