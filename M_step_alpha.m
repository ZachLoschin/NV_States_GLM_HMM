function [Pi_new, A_new, E_new, weightDiff] = M_step_alpha(input, output, eta_11, eta_12, eta_21, eta_22, gamma1, gamma2, E)
    % Calculate the pi_k initial state probability updates
    % Pi_k = gamma(Z_1k) / Sum{j}[gamma(Z1j)]
    Pi_new = [];
    
    % Calculate the new initiation probabilities
    sum_term = gamma1(1) + gamma2(1);
    Pi_new(1) = gamma1(1) / sum_term;
    Pi_new(2) = gamma2(1) / sum_term;
    
    
    % Calculate the new transition matrix using the eta values
    % # Expected transitions j->k / # expected transitions j->any state
    A_new(1,1) = sum(eta_11) / (sum(eta_11) + sum(eta_12));
    A_new(1,2) = sum(eta_12) / (sum(eta_11) + sum(eta_12));
    A_new(2,1) = sum(eta_21) / (sum(eta_21) + sum(eta_22));
    A_new(2,2) = sum(eta_22) / (sum(eta_21) + sum(eta_22));
    
    % -- Fit the new IRFs
    
    % Calcualte the true IRFs for each time window
    
    % Calculate all of the alpha function predictions
    [examples, window] = size(input);
    predicted_alphas = zeros(examples, window);

    for chim=1:length(input)
        predicted_alphas(chim, :) = ifft(fft(output(chim,:)) ./ fft(input(chim,:)));
    end

    

    disp("hi")
    % Normalize the weights so each timepoint sums to one across the
    % class responsabilities

    % --- WEIGHTING CORRECTLY --- %

    % Normalize the weights
    sum_gamma_term = ((gamma1) + (gamma2));
    weights1 = (gamma1) ./ sum_gamma_term;
    weights2 = (gamma2) ./ sum_gamma_term;

    % Track diff in weights, should diverge!
    weightDiff = mean(abs(weights1 - weights2));

    % Calculate new alpha functions
    % Combine predicted alpha functions for class 1
    class1_alpha = (predicted_alphas' * weights1) ./ sum(weights1);

    % Combine predicted alpha functions for class 2
    class2_alpha = (predicted_alphas' * weights2) ./ sum(weights2);

    E_new = [class1_alpha'; class2_alpha'];


    
    
    % % --- VITERBI WEIGHTING --- %
    % 
    % 
    % sum_gamma_term = ((gamma1) + (gamma2));
    % weights1 = (gamma1) ./ sum_gamma_term;
    % weights2 = (gamma2) ./ sum_gamma_term;
    % % 
    % weightDiff = mean(abs(weights1 - weights2));
    % % 
    % w1 = zeros(length(weights1), 1);
    % w2 = zeros(length(weights2), 1);
    % 
    % w1_ind = find(weights1>=weights2);
    % w2_ind = find(weights2>=weights1);
    % 
    % w1(w1_ind) = 1;
    % w2(w2_ind) = 1;
    % 
    % % Combine predicted alpha functions for class 1
    % class1_alpha = sum(predicted_alphas .* w1, 1) ./ length(w1_ind);
    % 
    % % Combine predicted alpha functions for class 2
    % class2_alpha = sum(predicted_alphas .* w2, 1) ./ length(w2_ind);
    % 
    % E_new = [class1_alpha; class2_alpha];
end

