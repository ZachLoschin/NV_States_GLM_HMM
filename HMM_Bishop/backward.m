function [B_norm] = backward(HMM, X, Cn)
    % Function to calculate the backward probabilities for an IOHMM
    % using the same scaling factors as the forward algorithm

    % Initialize variables
    [n, del] = size(X);  % Length of chain
    k = length(HMM.A(1,:));  % Number of states

    % Initialize backward probabilities
    B_norm = zeros(n, k);

    % Initialize backward probabilities at the last time step based on emission probabilities
    for state = 1:k
        B_norm(n, state) = 0.5;
    end

    % Iterate backward from second-to-last time step to the first time step
    for t = n-1:-1:1  % For each time point from second-to-last to the first
        for state = 1:k  % For each state 1->k
            B_current = zeros(1, k);

            % Compute backward probabilities for each next state
            for next_state = 1:k
                B_current(next_state) = B_norm(t+1, next_state) * HMM.A(state, next_state) * mvnpdf(X(t+1,:), HMM.U{next_state}, HMM.Sigma{next_state});
            end

            % Normalize backward probabilities using the same scaling factors as the forward algorithm
            B_norm(t, state) = sum(B_current) ./ Cn(t+1);
        end
    end
end



