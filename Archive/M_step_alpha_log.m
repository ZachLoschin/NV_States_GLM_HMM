function [Pi_new, A_new, E_new, weightDiff, state_prediction, weights1, weights2] = M_step_alpha_log(input, output, xi, gamma, E)
    % Calculate the pi_k initial state probability updates
    % Pi_k = gamma(Z_1k) / Sum{j}[gamma(Z1j)]
    Pi_new = [];
    
    % Calculate the new initiation probabilities
    denom = sum(exp(gamma(1,:)));
    Pi_new = exp(gamma(1,:)) ./ denom;
    
    % Calculate new transition probabilities
    
    % -- Fit the new IRFs
    [K, l] = size(E);
    N = length(input);

    A_new = zeros(2,2);
    for i = 1:K
        for j = 1:K
            A_new(i, j) = exp(log(sum(exp(xi(i,j,:)))) - log(sum(exp(gamma(1:N-1, i)))));
        end
    end

    % Normalize rows of the transition matrix
    sum_term = sum(A_new, 2);

    A_new(1,:) = A_new(1,:) ./ sum_term(1);
    A_new(2,:) = A_new(2,:) ./ sum_term(2);

    
    % UPDATE EMISSION IRFS %

    % Calculate all of the alpha function predictions
    [examples, window] = size(input);
    predicted_alphas = zeros(examples, window);

    [N, del] = size(input);

    for chim=1:N
        predicted_alphas(chim, :) = ifft(fft(output(chim,:)) ./ fft(input(chim,:)));
    end

    % Normalize the weights so each timepoint sums to one across the
    % class responsabilities

    % --- WEIGHTING CORRECTLY --- %

    % Normalize the weights
    weights = exp(gamma);
    sum_term = sum(weights, 2);
    weights = weights ./ sum_term;

    % Track diff in weights, should diverge!
    weightDiff = mean(abs(weights(:,1) - weights(:,2)));

    % Calculate new alpha functions
    % Combine predicted alpha functions for class 1
    class1_alpha = (predicted_alphas' * weights(:,1)) ./ sum(weights(:,1));

    % Combine predicted alpha functions for class 2
    class2_alpha = (predicted_alphas' * weights(:,2)) ./ sum(weights(:,2));

    E_new = [class1_alpha'; class2_alpha'];

    % -- Calculate each states prediction based on max weight
    w1 = zeros(length(weights(:,1)), 1) + 2;
    w1_ind = find(weights(:,1)>=weights(:,2));
    w1(w1_ind) = 1;

    weights1 = weights(:,1);
    weights2 = weights(:,2);

    state_prediction = w1;
end

