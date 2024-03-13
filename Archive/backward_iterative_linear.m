function logB = backward_iterative_linear(x, u, A, E)
    % Function that calculates the backward log probabilities for an IOHMM

   % Inputs:   x = Nx70 vector of observed outputs
    %           u = Nx70 vector of inputs
    %           A = KxK matrix of transition probabilities
    %           E = 2xK matrix of emission probabilities (slope, intercept)
    
    % Initialize variables
    [n, del] = size(x);  % Length of chain
    k = length(A(1,:));  % Number of states

    % Starting points of iterations
    logB = zeros(n, k);  % Initialize log probabilities matrix

    % Initialize backward probabilities at the last time step based on emission probabilities
    for state = 1:k
        logB(n, state) = 0;  % initialize to 0 -> 1 in normal space
    end

    % Iterate backward from second-to-last time step to the first time step
    for t = n-1:-1:1  % For each time point from second-to-last to the first
        for state = 1:k  % For each state 1->k
            values_to_sum = zeros(1, k);
            for next_state = 1:k
                values_to_sum(next_state) = logB(t+1, next_state) + log(A(state, next_state)) + linear_emit(x(t+1, :), u(t+1, :), E(next_state, :));
            end
            logB(t, state) = logsumexp(values_to_sum); % Calculate the sum of exponentials using logsumexp
        end
    end
end



