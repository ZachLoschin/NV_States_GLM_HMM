% Zachary Loschinskey
% Anna Devor Rotation
% January 2024
% Neurovascular States Project

% Function to calculate the forward probability of an IOHMM with dual alpha
% function defined relationships between the input and output.

function logF_lse = forward_iterative_linear(x, u, A, E, pi)
    % Function that calculates the forward log probabilities for an IOHMM

    % Inputs:   x = Nx70 vector of observed outputs
    %           u = Nx70 vector of inputs
    %           A = KxK matrix of transition probabilities
    %           E = 2xK matrix of emission probabilities

    % Initialize variables
    [n, del] = size(x);  % Length of chain I SWEAR
    k = length(A(1,:));  % Number of states

    % Starting points of iterations
    logF_lse = zeros(n, k);  % Initialize log probabilities matrix

    for state = 1:k
        logF_lse(1, state) = log(pi(k)) + linear_emit(x(1, :), u(1, :), E(state,:));
    end

    for t = 2:n  % For each time point 2->n
        for state = 1:k  % For each state 1->k
            values_to_sum = zeros(1, k);
            for prev_state = 1:k
                values_to_sum(prev_state) = logF_lse(t-1, prev_state) + log(A(prev_state, state));
            end
            logF_lse(t, state) = logsumexp(values_to_sum)  + linear_emit(x(t, :), u(t, :), E(state, :)); % Logsumexp to avoid underflow
        end
    end
end
