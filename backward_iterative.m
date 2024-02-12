function [B_norm] = backward_iterative(x, u, A, E, Cn)
    % Function to calculate the backward probabilities for an IOHMM
    % using the same scaling factors as the forward algorithm

    % Inputs:   x = Nx70 vector of observed outputs
    %           u = Nx70 vector of inputs
    %           A = KxK matrix of transition probabilities
    %           E = 2xK matrix of emission probabilities
    %           Cn = 1xn vector of scaling factors from the forward algorithm

    % Initialize variables
    n = length(x);  % Length of chain
    k = length(A(1,:));  % Number of states

    % Initialize backward probabilities
    B_norm = zeros(n, k);

    % Initialize backward probabilities at the last time step based on emission probabilities
    for state = 1:k
        B_norm(n, state) = 0.5; % Initialize to 0.5
    end

    % Iterate backward from second-to-last time step to the first time step
    for t = n-1:-1:1  % For each time point from second-to-last to the first
        for state = 1:k  % For each state 1->k
            B_current = zeros(1, k);

            % Compute backward probabilities for each next state
            for next_state = 1:k
                B_current(next_state) = B_norm(t+1, next_state) * A(state, next_state) * alpha_emit_prob(x(t+1, :), u(t+1, :), E(next_state, :));
            end

            % Normalize backward probabilities using the same scaling factors as the forward algorithm
            B_norm(t, state) = sum(B_current) ./ Cn(t+1);
        end
    end
end



