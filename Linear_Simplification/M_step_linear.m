function [Pi_new, A_new, E_new, gamma1, gamma2, weights1, weights2] = M_step_linear(input, output, eta_11, eta_12, eta_21, eta_22, gamma1, gamma2, E)
    % Calculate the pi_k initial state probability updates
    % Pi_k = gamma(Z_1k) / Sum{j}[gamma(Z1j)]
    Pi_new = [];
    
    % Calculate the new initiation probabilities
    sum_term = gamma1(1) + gamma2(1);
    Pi_new(1) = gamma1(1) / sum_term;
    Pi_new(2) = gamma2(1) / sum_term;
    
    
    % Calculate the new transition matrix using the eta values
    % # Expected transitions j->k / # expected transitions j->any state
    A_new(1,1) = sum(eta_11(2:end)) / (sum(eta_11(2:end)) + sum(eta_12(2:end)));
    A_new(1,2) = sum(eta_12(2:end)) / (sum(eta_11(2:end)) + sum(eta_12(2:end)));
    A_new(2,1) = sum(eta_21(2:end)) / (sum(eta_21(2:end)) + sum(eta_22(2:end)));
    A_new(2,2) = sum(eta_22(2:end)) / (sum(eta_21(2:end)) + sum(eta_22(2:end)));
    
    % -- Fit the new IRFs
    
    % Define the linear model
    linearModel = fittype('a*x + b', 'coefficients', {'a', 'b'});
    
    weights1 = gamma1;
    weights2 = gamma2;
    
    % -- MAXES
    % Find the indices where weights1 is greater than weights2
    greater_indices = weights1 > weights2;
    lesser_indices = weights1 < weights2;
    
    % Create a vector of zeros with the same length as weights1 and weights2
    w1 = zeros(size(weights1));
    
    % Set the elements at the indices where weights1 > weights2 to 1
    w1(greater_indices) = 1;
    w2 = ~w1;

    w1 = repelem(w1, 70);
    w2 = double(repelem(w2, 70));

    weights1_fit = repelem(weights1, 70);
    weights2_fit = repelem(weights2, 70);

    % Give the current guess as the fit parameters
    startPoint1 = E(1, :);
    startPoint2 = E(2, :);

    % Transform the data into correct shape for fit()
    input_fit = input(:);
    output_fit = output(:);

    % Repeat the element weights to make same size
    weights1_fit = repelem(weights1, 70);
    weights2_fit = repelem(weights2, 70);

    % Create a fit object with starting points
    FO1 = fit(input_fit(greater_indices), output_fit(greater_indices), linearModel);
    FO2 = fit(input_fit(lesser_indices), output_fit(lesser_indices), linearModel);

    % Exctract the parameters from the fit object
    E_new = [FO1.a FO1.b; FO2.a FO2.b];

    % -- visualization of weighting
    input_1 = input(greater_indices);
    output_1 = output(greater_indices);

    input_2 = input(lesser_indices);
    output_2 = output(lesser_indices);
    % 
    % figure()
    % scatter(input_1, output_1, 'r')
    % hold on
    % scatter(input_2, output_2, 'b')
    % legend("Class 1 Max, Class 2 Max")
    % xlabel("x")
    % ylabel("y")


end

